use bstr::ByteSlice;

pub fn extract_name_bytes<'a>(id: &'a str, sep: &Option<String>) -> &'a [u8] {
    let id = id.as_bytes();
    if let Some(sep) = sep {
        if let Some(pos) = id.find(sep) {
            return &id[..pos];
        }
    }
    id
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_name_bytes() {
        assert_eq!(
            extract_name_bytes("foo|bar|baz", &Some("|".to_owned())),
            b"foo"
        );
        assert_eq!(
            extract_name_bytes("fooabbarabbaz", &Some("ab".to_owned())),
            b"foo"
        );
    }
}
