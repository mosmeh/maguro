mod rank9b;
mod suffix_array;

use crate::{sequence, utils};
use bio::io::fasta::{self, FastaRead};
use bitvec::prelude::*;
use rank9b::Rank9b;
use serde::{Deserialize, Serialize};
use std::{io, path::Path};
use suffix_array::SuffixArray;

pub const DELIMITER: u8 = b'$';

#[derive(Serialize, Deserialize)]
pub struct Index {
    pub(crate) seq: Vec<u8>,
    ends: Vec<usize>,
    rank_dict: Rank9b,
    name_arena: Vec<u8>,
    name_ends: Vec<usize>,
    pub(crate) sa: SuffixArray,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SequenceId(pub usize);

impl Index {
    pub fn num_seqs(&self) -> usize {
        self.ends.len() - 1
    }

    pub fn seq_name(&self, seq_id: SequenceId) -> &[u8] {
        &self.name_arena[self.name_ends[seq_id.0]..self.name_ends[seq_id.0 + 1] - 1]
    }

    pub fn seq(&self, seq_id: SequenceId) -> &[u8] {
        &self.seq[self.seq_range(seq_id)]
    }

    pub(crate) fn seq_id_from_pos(&self, pos: usize) -> SequenceId {
        assert!(self.ends[0] <= pos && pos < self.seq.len(), "Out of bounds");

        let rank = self.rank_dict.rank(pos) as usize;
        assert!(rank >= 1);

        let seq_id = rank - 1;
        assert!(seq_id < self.num_seqs());

        SequenceId(seq_id)
    }

    pub(crate) fn seq_range(&self, seq_id: SequenceId) -> std::ops::Range<usize> {
        self.ends[seq_id.0]..(self.ends[seq_id.0 + 1] - 1)
    }
}

pub struct IndexBuilder<R: io::Read> {
    reader: fasta::Reader<R>,
    bucket_width: usize,
    header_sep: Option<String>,
}

impl IndexBuilder<std::fs::File> {
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self::from_reader(fasta::Reader::from_file(path)?))
    }
}

impl<R: io::Read> IndexBuilder<R> {
    fn from_reader(reader: fasta::Reader<R>) -> Self {
        Self {
            reader,
            bucket_width: 8,
            header_sep: None,
        }
    }

    pub fn new(reader: R) -> Self {
        Self::from_reader(fasta::Reader::new(reader))
    }

    pub fn bucket_width(mut self, bucket_width: usize) -> Self {
        self.bucket_width = bucket_width;
        self
    }

    pub fn header_sep(mut self, header_sep: String) -> Self {
        self.header_sep = Some(header_sep);
        self
    }

    pub fn build(mut self) -> io::Result<Index> {
        let mut seq = vec![DELIMITER];
        let mut ends = vec![1];
        let mut name_arena = Vec::new();
        let mut name_ends = vec![0];

        let mut record = fasta::Record::new();
        self.reader.read(&mut record)?;
        while !record.is_empty() {
            record.check().unwrap();

            seq.extend_from_slice(record.seq());
            seq.push(DELIMITER);
            ends.push(seq.len());

            name_arena.extend(utils::extract_byte_name(record.id(), &self.header_sep));
            name_arena.push(b'\n');
            name_ends.push(name_arena.len());

            self.reader.read(&mut record)?;
        }

        sequence::encode_in_place(&mut seq);

        let mut bvec: BitVec<Lsb0, u64> = BitVec::new();
        bvec.resize(seq.len(), false);
        for end in &ends {
            *bvec.get_mut(end - 1).unwrap() = true;
        }

        let sa = SuffixArray::new(&seq, self.bucket_width);
        Ok(Index {
            seq,
            ends,
            rank_dict: Rank9b::from_bit_vec(bvec),
            name_arena,
            name_ends,
            sa,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_seq_id(seqs: &[&[u8]], pos: &[usize], expected: &[usize]) {
        let mut fasta: Vec<u8> = Vec::new();
        for seq in seqs {
            fasta.extend(b">foo\n");
            fasta.extend(*seq);
            fasta.push(b'\n');
        }

        let cursor = std::io::Cursor::new(&fasta);
        let index = IndexBuilder::new(cursor).build().unwrap();

        let got: Vec<_> = pos
            .iter()
            .map(|pos| index.seq_id_from_pos(*pos).0)
            .collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn seq_id() {
        check_seq_id(&[b"agctagt"], &[1, 3, 5, 7, 8], &[0; 5]);
        check_seq_id(
            &[b"agct", b"tgta"],
            &[1, 4, 5, 6, 9, 10],
            &[0, 0, 0, 1, 1, 1],
        );
        check_seq_id(
            &[b"atcgggatatatggagagcttagag", b"tttagagggttcttcgggatt"],
            &[1, 10, 25, 26, 27, 35, 47, 48],
            &[0, 0, 0, 0, 1, 1, 1, 1],
        );
    }

    #[test]
    #[should_panic(expected = "Out of bounds")]
    fn out_of_bounds_left() {
        check_seq_id(&[b"agctagt"], &[0], &[0]);
    }

    #[test]
    #[should_panic(expected = "Out of bounds")]
    fn out_of_bounds_right() {
        check_seq_id(&[b"agctagt"], &[9], &[0]);
    }
}
