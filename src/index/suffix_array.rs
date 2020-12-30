use crate::sequence;
use bio::data_structures::suffix_array::suffix_array;
use serde::{Deserialize, Serialize};
use std::ops::Range;

#[derive(Serialize, Deserialize)]
pub struct SuffixArray {
    pub array: Vec<usize>,
    child: Vec<usize>,
    buckets: Vec<Range<usize>>,
    bucket_width: usize,
}

impl SuffixArray {
    pub fn new(text: &[u8], bucket_width: usize) -> Self {
        assert!(bucket_width * 2 < std::mem::size_of::<usize>() * 8);

        let array = suffix_array(text);
        let lcp = make_lcp_array(text, &array);
        let child = make_child_table(&lcp);

        let bucket_len = 1 << (2 * bucket_width);
        let mut prefix = vec![0; bucket_width];
        let mut buckets = Vec::with_capacity(bucket_len);

        for i in 0..bucket_len {
            for (bit, x) in prefix.iter_mut().enumerate() {
                *x = sequence::two_bit_to_code((i >> (2 * bit)) as u8);
            }
            buckets.push(
                search_suffix_array(&array, &child, text, &prefix, 0, 0, text.len(), 0)
                    .unwrap_or_default(),
            );
        }

        Self {
            array,
            child,
            buckets,
            bucket_width,
        }
    }

    pub fn search(&self, text: &[u8], query: &[u8]) -> Option<Range<usize>> {
        if query.len() < self.bucket_width {
            return search_suffix_array(&self.array, &self.child, text, query, 0, 0, text.len(), 0);
        }

        let mut idx = 0;
        for (i, x) in query[..self.bucket_width].iter().enumerate() {
            idx |= (sequence::code_to_two_bit(*x) as usize) << (2 * i);
        }

        let range = &self.buckets[idx];
        if range.start == range.end {
            return None;
        } else if query.len() == self.bucket_width {
            return Some(range.clone());
        }

        let store_pos = if self.child[range.start] < range.end {
            range.start
        } else {
            range.end - 1
        };

        search_suffix_array(
            &self.array,
            &self.child,
            text,
            query,
            self.bucket_width,
            range.start,
            range.end,
            store_pos,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn search_suffix_array(
    array: &[usize],
    child: &[usize],
    text: &[u8],
    query: &[u8],
    mut depth: usize,
    mut begin: usize,
    mut end: usize,
    mut store_pos: usize,
) -> Option<Range<usize>> {
    'label0: loop {
        if depth == query.len() {
            return Some(begin..end);
        }
        let q = query[depth];
        let mut b = text[array[begin] + depth];
        if q < b {
            return None;
        }

        'label1: loop {
            let e = text[array[end - 1] + depth];
            if q > e {
                return None;
            }

            loop {
                if b == e {
                    depth += 1;
                    continue 'label0;
                }
                let mid = child[store_pos];
                let m = text[array[mid] + depth];
                if q < m {
                    end = mid;
                    store_pos = mid - 1;
                    continue 'label1;
                }
                begin = mid;
                b = m;
                store_pos = mid;
            }
        }
    }
}

fn make_lcp_array(text: &[u8], sa: &[usize]) -> Vec<usize> {
    let len = text.len();
    let mut inv_sa = vec![0; len];
    let mut lcp = vec![0; len];
    for i in 0..len {
        inv_sa[sa[i]] = i;
    }

    // Kasai's algorithm
    let mut k = 0;
    for i in 0..len {
        if inv_sa[i] == len - 1 {
            k = 0;
            continue;
        }
        let j = sa[inv_sa[i] + 1];
        while i.max(j) + k < len && text[i + k] == text[j + k] {
            k += 1;
        }
        lcp[inv_sa[i] + 1] = k;
        if k > 0 {
            k -= 1;
        }
    }

    lcp
}

fn make_child_table(lcp: &[usize]) -> Vec<usize> {
    let mut child_table = vec![0; lcp.len()];
    let mut stack = vec![(0, lcp.len(), 0)];
    while let Some((beg, end, store_pos)) = stack.pop() {
        if end - beg < 2 {
            continue;
        }
        let (mid, _) = lcp
            .iter()
            .enumerate()
            .skip(beg + 1)
            .take(end - beg - 1)
            .min_by(|(_, x), (_, y)| x.cmp(&y))
            .unwrap();
        child_table[store_pos] = mid;
        stack.push((beg, mid, mid - 1));
        stack.push((mid, end, mid))
    }
    child_table
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence;
    use itertools::Itertools;

    fn check_search(text: &[u8], query: &[u8], expected: Option<&[usize]>) {
        let text = sequence::encode(text);
        let query = sequence::encode(query);
        let sa = SuffixArray::new(&text, 2);
        let got: Option<Vec<usize>> = sa
            .search(&text, &query)
            .map(|result| result.map(|i| sa.array[i]).sorted().collect());
        if expected.is_some() {
            assert_eq!(got.unwrap(), expected.unwrap());
        } else {
            assert!(got.is_none());
        }
    }

    #[test]
    fn search() {
        check_search(b"agctagt$", b"gct", Some(&[1]));
        check_search(b"agctagt$", b"agctagta", None);
        check_search(b"agctagt$", b"gca", None);
        check_search(b"atcgggatatatggagagcttagag$", b"gag", Some(&[13, 15, 22]));
        check_search(b"agct$tgta$", b"a", Some(&[0, 8]));
        check_search(
            b"atcgggatatatggagagcttagag$tttagagggttcttcgggatt$",
            b"gag",
            Some(&[13, 15, 22, 30]),
        );
    }
}
