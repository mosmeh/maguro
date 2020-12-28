use bstr::ByteSlice;
use std::ops::Range;

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
