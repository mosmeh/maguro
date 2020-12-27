use bio::data_structures::suffix_array::suffix_array;
use serde::{Deserialize, Serialize};
use std::ops::Range;

#[derive(Serialize, Deserialize)]
pub struct SuffixArray {
    pub array: Vec<usize>,
    child: Vec<usize>,
}

impl SuffixArray {
    pub fn new(text: &[u8]) -> Self {
        let array = suffix_array(text);
        let lcp = make_lcp_array(text, &array);
        Self {
            array,
            child: make_child_table(&lcp),
        }
    }

    pub fn search(&self, text: &[u8], query: &[u8]) -> Option<Range<usize>> {
        let mut depth = 0;
        let mut begin = 0;
        let mut end = text.len();
        let mut store_pos = 0;
        'label0: loop {
            if depth == query.len() {
                return Some(begin..end);
            }
            let q = query[depth];
            let mut b = text[self.array[begin] + depth];
            if q < b {
                return None;
            }

            'label1: loop {
                let e = text[self.array[end - 1] + depth];
                if q > e {
                    return None;
                }

                loop {
                    if b == e {
                        depth += 1;
                        continue 'label0;
                    }
                    let mid = self.child[store_pos];
                    let m = text[self.array[mid] + depth];
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
    use itertools::Itertools;

    fn check_search(text: &[u8], query: &[u8], expected: Option<&[usize]>) {
        let sa = SuffixArray::new(text);
        let got: Option<Vec<usize>> = sa
            .search(text, query)
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
