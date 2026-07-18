//! Input handling: transparent gzip/bgzip decompression and stdin, mapped to memory.
//!
//! The two-pass architecture needs random access, which a pipe or a compressed stream can't
//! provide, so compressed or piped input is spooled to a temp file that is then memory-mapped
//! and removed on drop.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read};
use std::path::PathBuf;

use flate2::read::MultiGzDecoder;
use memmap2::Mmap;

/// A memory-mapped input, owning any temp file created for decompressed/piped data.
pub struct MappedInput {
    pub mmap: Mmap,
    temp: Option<PathBuf>,
}

impl Drop for MappedInput {
    fn drop(&mut self) {
        if let Some(p) = &self.temp {
            let _ = std::fs::remove_file(p);
        }
    }
}

fn temp_path() -> PathBuf {
    std::env::temp_dir().join(format!("snpick-{}.tmp", std::process::id()))
}

fn spool_to_temp(mut reader: impl Read) -> io::Result<(Mmap, PathBuf)> {
    let tp = temp_path();
    {
        let mut out = BufWriter::new(File::create(&tp)?);
        io::copy(&mut reader, &mut out)?;
    }
    let f = File::open(&tp)?;
    let mmap = unsafe { Mmap::map(&f)? };
    Ok((mmap, tp))
}

/// Map the input for reading. `path == "-"` reads stdin; gzip/bgzip input (detected by magic
/// bytes) is decompressed transparently.
pub fn map_input(path: &str) -> io::Result<MappedInput> {
    if path == "-" {
        // Peek the first two bytes, then chain them back so gzip over stdin is detected.
        let mut reader = BufReader::new(io::stdin().lock());
        let mut magic = [0u8; 2];
        let n = reader.read(&mut magic)?;
        let head = std::io::Cursor::new(magic[..n].to_vec());
        let stream = head.chain(reader);
        let (mmap, tp) = if n == 2 && magic == [0x1f, 0x8b] {
            spool_to_temp(MultiGzDecoder::new(stream))?
        } else {
            spool_to_temp(stream)?
        };
        return Ok(MappedInput { mmap, temp: Some(tp) });
    }

    let f = File::open(path)
        .map_err(|e| io::Error::new(e.kind(), format!("Cannot open '{}': {}", path, e)))?;

    // Peek the first two bytes for the gzip magic (0x1f 0x8b); bgzip is a valid gzip stream.
    let mut magic = [0u8; 2];
    let n = (&f).read(&mut magic)?;
    let is_gzip = n == 2 && magic == [0x1f, 0x8b];

    if is_gzip {
        let inf = File::open(path)?;
        let (mmap, tp) = spool_to_temp(MultiGzDecoder::new(BufReader::new(inf)))?;
        Ok(MappedInput { mmap, temp: Some(tp) })
    } else {
        let mmap = unsafe { Mmap::map(&f)? };
        Ok(MappedInput { mmap, temp: None })
    }
}
