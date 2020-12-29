use bstr::ByteSlice;

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
