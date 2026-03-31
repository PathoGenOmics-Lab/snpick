//! Zero-copy FASTA parser over memory-mapped data.
//!
//! Records are indexed by scanning for `>` headers and tracking sequence offsets.
//! No data is copied — IDs and descriptions are `&[u8]` slices into the mmap.

use std::io;

use crate::types::{SeqLayout, MAX_SEQ_LENGTH};

/// A FASTA record as zero-copy slices into memory-mapped data.
pub struct FastaRecord<'a> {
    pub id: &'a [u8],
    pub desc: &'a [u8],
    pub seq_offset: usize,
}

/// Index FASTA records from memory-mapped data. Zero-copy: stores `&[u8]` slices.
///
/// Returns `(records, seq_length, layout)`.
/// All sequences must have the same length (alignment requirement).
pub fn index_fasta(data: &[u8]) -> io::Result<(Vec<FastaRecord<'_>>, usize, SeqLayout)> {
    if data.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Input FASTA is empty (0 bytes)."));
    }

    let mut records = Vec::new();
    let mut pos = 0;
    let len = data.len();
    let mut seq_length = 0usize;
    let mut is_single_line = true;

    while pos < len {
        while pos < len && (data[pos] == b'\n' || data[pos] == b'\r') { pos += 1; }
        if pos >= len { break; }

        if data[pos] != b'>' {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Expected '>' at byte {}, got '{}'.", pos, data[pos] as char)));
        }
        pos += 1;

        // Header line
        let header_start = pos;
        while pos < len && data[pos] != b'\n' && data[pos] != b'\r' { pos += 1; }
        let header = &data[header_start..pos];
        if pos < len && data[pos] == b'\r' { pos += 1; }
        if pos < len && data[pos] == b'\n' { pos += 1; }

        let (id, desc) = if let Some(sp) = header.iter().position(|&b| b == b' ' || b == b'\t') {
            (&header[..sp], &header[sp + 1..])
        } else {
            (header, &data[0..0])
        };

        let seq_offset = pos;

        // Count bases
        let mut seq_len = 0usize;
        let mut line_count = 0usize;
        while pos < len && data[pos] != b'>' {
            let line_start = pos;
            while pos < len && data[pos] != b'\n' && data[pos] != b'\r' { pos += 1; }
            let line_len = pos - line_start;
            if line_len > 0 {
                seq_len += line_len;
                line_count += 1;
            }
            if pos < len && data[pos] == b'\r' { pos += 1; }
            if pos < len && data[pos] == b'\n' { pos += 1; }
        }

        if line_count > 1 { is_single_line = false; }

        if records.is_empty() {
            if seq_len == 0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "First sequence has length 0."));
            }
            if seq_len > MAX_SEQ_LENGTH {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    format!("Sequence length {} exceeds maximum.", seq_len)));
            }
            seq_length = seq_len;
        } else if seq_len != seq_length {
            let id_str = std::str::from_utf8(id).unwrap_or("?");
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Sequence '{}' (#{}) has length {} but expected {}.",
                    id_str, records.len() + 1, seq_len, seq_length)));
        }

        records.push(FastaRecord { id, desc, seq_offset });
    }

    if records.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Input FASTA is empty."));
    }

    Ok((records, seq_length, SeqLayout { single_line: is_single_line }))
}

/// Extract reference sequence from first record.
pub fn get_ref_seq(data: &[u8], rec: &FastaRecord, seq_length: usize, layout: SeqLayout) -> Vec<u8> {
    if layout.single_line {
        data[rec.seq_offset..rec.seq_offset + seq_length].to_vec()
    } else {
        let mut seq = Vec::with_capacity(seq_length);
        let mut pos = rec.seq_offset;
        let end = data.len();
        while seq.len() < seq_length && pos < end {
            let b = data[pos];
            pos += 1;
            if b != b'\n' && b != b'\r' { seq.push(b); }
        }
        seq
    }
}
