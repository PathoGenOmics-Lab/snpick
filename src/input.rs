//! Input handling: transparent gzip/bgzip decompression and stdin, mapped to memory.
//!
//! The two-pass architecture needs random access, which a pipe or a compressed stream can't
//! provide, so compressed or piped input is spooled to a temp file that is then memory-mapped
//! and removed on drop.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
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

/// Read up to 2 magic bytes, tolerating short reads (a pipe's first read may return < 2 bytes).
fn read_magic(mut r: impl Read) -> io::Result<(usize, [u8; 2])> {
    let mut buf = [0u8; 2];
    let mut n = 0;
    while n < 2 {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    Ok((n, buf))
}

fn spool_to_temp(mut reader: impl Read) -> io::Result<(Mmap, PathBuf)> {
    let tp = temp_path();
    let build = (|| -> io::Result<Mmap> {
        {
            let mut out = BufWriter::new(File::create(&tp)?);
            io::copy(&mut reader, &mut out)?;
            out.flush()?;
        }
        let f = File::open(&tp)?;
        unsafe { Mmap::map(&f) }
    })();
    match build {
        Ok(mmap) => Ok((mmap, tp)),
        Err(e) => {
            // Don't leak the spool file if decompression / mmap failed.
            let _ = std::fs::remove_file(&tp);
            Err(e)
        }
    }
}

/// Map the input for reading. `path == "-"` reads stdin; gzip/bgzip input (detected by magic
/// bytes) is decompressed transparently.
pub fn map_input(path: &str) -> io::Result<MappedInput> {
    if path == "-" {
        // Peek the first two bytes, then chain them back so gzip over stdin is detected.
        let mut reader = BufReader::new(io::stdin().lock());
        let (n, magic) = read_magic(&mut reader)?;
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
    let (n, magic) = read_magic(&f)?;
    let is_gzip = n == 2 && magic == [0x1f, 0x8b];

    if is_gzip {
        // Chain the peeked magic back onto the same handle rather than reopening the path.
        // Reopening rewinds a regular file but NOT a non-seekable path (process substitution
        // `<(gzip -c ...)` or a FIFO), where the magic bytes are already consumed — the decoder
        // would then start mid-stream and fail with "invalid gzip header". This mirrors the
        // stdin branch and works for both seekable and non-seekable inputs.
        let head = std::io::Cursor::new(magic[..n].to_vec());
        let stream = head.chain(BufReader::new(f));
        let (mmap, tp) = spool_to_temp(MultiGzDecoder::new(stream))?;
        Ok(MappedInput { mmap, temp: Some(tp) })
    } else {
        // Mmap maps from offset 0 regardless of the read position left by read_magic.
        let mmap = unsafe { Mmap::map(&f)? };
        Ok(MappedInput { mmap, temp: None })
    }
}

#[cfg(test)]
mod tests {
    use super::map_input;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs::File;
    use std::io::Write;

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut e = GzEncoder::new(Vec::new(), Compression::fast());
        e.write_all(bytes).unwrap();
        e.finish().unwrap()
    }

    #[test]
    fn gzip_regular_file_decodes() {
        let raw = b">s1\nATGC\n>s2\nATCC\n";
        let p = std::env::temp_dir().join(format!("snpick-gz-{}.fa.gz", std::process::id()));
        std::fs::write(&p, gzip(raw)).unwrap();
        let m = map_input(p.to_str().unwrap()).unwrap();
        assert_eq!(&m.mmap[..], raw);
        std::fs::remove_file(&p).ok();
    }

    // A gzip stream delivered through a non-seekable FIFO must still decode. On the old code the
    // gzip branch reopened the path, which on a pipe yields a fresh handle whose magic bytes were
    // already consumed, so decoding failed with "invalid gzip header".
    #[cfg(unix)]
    #[test]
    fn gzip_over_fifo_is_rewind_safe() {
        use std::process::Command;
        let raw = b">s1\nATGC\n>s2\nATCC\n";
        let gz = gzip(raw);
        let fifo = std::env::temp_dir().join(format!("snpick-fifo-{}.gz", std::process::id()));
        let _ = std::fs::remove_file(&fifo);
        assert!(Command::new("mkfifo").arg(&fifo).status().unwrap().success());

        let wpath = fifo.clone();
        let writer = std::thread::spawn(move || {
            // File::create on a FIFO blocks until the reader opens; then stream and close.
            if let Ok(mut f) = File::create(&wpath) {
                let _ = f.write_all(&gz);
            }
        });
        let m = map_input(fifo.to_str().unwrap()).unwrap();
        assert_eq!(&m.mmap[..], raw);
        writer.join().unwrap();
        std::fs::remove_file(&fifo).ok();
    }
}
