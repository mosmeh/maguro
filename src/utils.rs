use bstr::ByteSlice;
use std::ops::Range;

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

pub fn encode_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|x| ENCODE_TABLE[*x as usize]).collect()
}

pub fn encode_seq_in_place(seq: &mut [u8]) {
    for x in seq.iter_mut() {
        *x = ENCODE_TABLE[*x as usize];
    }
}

pub fn decode_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|x| DECODE_TABLE[*x as usize]).collect()
}

pub fn decode_seq_in_place(seq: &mut [u8]) {
    for x in seq.iter_mut() {
        *x = DECODE_TABLE[*x as usize];
    }
}

pub fn upper_bound<T: Ord>(slice: &[T], x: &T) -> usize {
    use std::cmp::Ordering::*;

    let mut left = 0;
    let mut right = slice.len();
    while left != right {
        let mid = left + (right - left) / 2;
        match slice[mid].cmp(x) {
            Less | Equal => {
                left = mid + 1;
            }
            Greater => {
                right = mid;
            }
        }
    }
    left
}

pub fn extract_byte_name<'a>(id: &'a str, sep: &Option<String>) -> &'a [u8] {
    let id = id.as_bytes();

    let split_pos = if let Some(sep) = sep {
        id.find(sep)
    } else {
        id.iter().position(u8::is_ascii_whitespace)
    };

    if let Some(i) = split_pos {
        &id[..i]
    } else {
        id
    }
}

pub fn intersect(a: &Range<usize>, b: &Range<usize>) -> Range<usize> {
    if b.start >= a.end || a.start >= b.end {
        return Default::default();
    }
    (a.start.max(b.start))..(a.end.max(b.end))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_byte_name() {
        assert_eq!(extract_byte_name("foo bar baz", &None), b"foo");
        assert_eq!(extract_byte_name("foo\tbar\tbaz", &None), b"foo");
        assert_eq!(
            extract_byte_name("foo|bar|baz", &Some("|".to_owned())),
            b"foo"
        );
        assert_eq!(
            extract_byte_name("fooabbarabbaz", &Some("ab".to_owned())),
            b"foo"
        );
    }
}
