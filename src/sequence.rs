pub const DUMMY_CODE: u8 = 5;

const ENCODE_TABLE: [u8; 256] = {
    let mut table = [DUMMY_CODE; 256];
    table[crate::index::DELIMITER as usize] = 0;
    table[b'A' as usize] = 1;
    table[b'a' as usize] = 1;
    table[b'C' as usize] = 2;
    table[b'c' as usize] = 2;
    table[b'G' as usize] = 3;
    table[b'g' as usize] = 3;
    table[b'T' as usize] = 4;
    table[b't' as usize] = 4;
    table
};

const DECODE_TABLE: [u8; 256] = {
    let mut table = [b'N'; 256];
    table[1] = b'A';
    table[2] = b'C';
    table[3] = b'G';
    table[4] = b'T';
    table
};

const COMPLEMENT_TABLE: [u8; 256] = {
    let mut table = [0; 256];
    table[1] = 4;
    table[2] = 3;
    table[3] = 2;
    table[4] = 1;
    table[5] = 5;
    table[b'A' as usize] = b'T';
    table[b'a' as usize] = b'T';
    table[b'C' as usize] = b'G';
    table[b'c' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'g' as usize] = b'C';
    table[b'T' as usize] = b'A';
    table[b't' as usize] = b'A';
    table
};

pub fn encode(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|x| ENCODE_TABLE[*x as usize]).collect()
}

pub fn encode_in_place(seq: &mut [u8]) {
    for x in seq.iter_mut() {
        *x = ENCODE_TABLE[*x as usize];
    }
}

pub fn decode(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|x| DECODE_TABLE[*x as usize]).collect()
}

pub fn decode_in_place(seq: &mut [u8]) {
    for x in seq.iter_mut() {
        *x = DECODE_TABLE[*x as usize];
    }
}

pub fn reverse(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().copied().collect()
}

pub fn complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|x| COMPLEMENT_TABLE[*x as usize]).collect()
}

pub fn complement_in_place(seq: &mut [u8]) {
    for x in seq.iter_mut() {
        *x = COMPLEMENT_TABLE[*x as usize];
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|x| COMPLEMENT_TABLE[*x as usize])
        .collect()
}

#[inline]
pub fn code_to_two_bit(x: u8) -> u8 {
    (x - 1) & 0b11
}

#[inline]
pub fn two_bit_to_code(x: u8) -> u8 {
    (x & 0b11) + 1
}
