use bitvec::prelude::*;
use serde::{Deserialize, Serialize};

// Adapted from https://github.com/foudrer/Sux

#[derive(Serialize, Deserialize)]
pub struct Rank9b {
    bits: Vec<u64>,
    counts: Vec<u64>,
}

impl Rank9b {
    pub fn from_vec(bits: Vec<u64>, num_bits: usize) -> Self {
        let num_words = (num_bits + 63) / 64;
        let num_counts = ((num_bits + 64 * 8 - 1) / (64 * 8)) * 2;

        let mut counts = vec![0; num_counts + 1];
        let mut c = 0u64;
        let mut pos = 0;
        for i in (0..num_words).step_by(8) {
            counts[pos] = c;
            c += bits[i].count_ones() as u64;
            for j in 1..8 {
                counts[pos + 1] |= (c - counts[pos]) << (63 - (9 * j));
                if i + j < num_words {
                    c += bits[i + j].count_ones() as u64;
                }
            }
            pos += 2;
        }

        counts[num_counts as usize] = c;
        assert!(c <= num_bits as u64);

        Self { bits, counts }
    }

    pub fn from_bit_vec(bits: BitVec<Lsb0, u64>) -> Self {
        let num_bits = bits.len();
        Self::from_vec(bits.into_vec(), num_bits)
    }

    pub fn bit(&self, k: usize) -> bool {
        self.bits[k / 64] & (1 << (k % 64)) != 0
    }

    pub fn rank(&self, k: usize) -> u64 {
        let word = k / 64;
        let block = (word / 4) & !1;
        let offset = word % 8;
        self.counts[block]
            + (self.counts[block + 1] >> (63 - offset * 9) & 0x1FF)
            + (self.bits[word] & ((1u64 << (k % 64)) - 1)).count_ones() as u64
    }

    pub fn size_bytes(&self) -> usize {
        std::mem::size_of::<u64>() * (self.bits.len() + self.counts.len())
    }
}
