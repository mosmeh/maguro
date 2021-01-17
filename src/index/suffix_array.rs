use crate::sequence;
use serde::{Deserialize, Serialize};
use std::ops::Range;
use sufsort_rs::sufsort::SA;

pub static STEP_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub static FOUND_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub static SEARCH_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

#[derive(Serialize, Deserialize)]
pub struct SuffixArray {
    pub array: Vec<u32>,
    pub child: Vec<u32>,
    pub buckets: Vec<Range<u32>>,
    pub bucket_width: usize,
}

impl SuffixArray {
    pub fn new(text: &[u8], bucket_width: usize) -> Self {
        assert!(text.len() <= u32::MAX as usize + 1);
        assert!(bucket_width * 2 < std::mem::size_of::<usize>() * 8);

        let sa = SA::<i32>::new(text);
        let array: Vec<_> = sa.sarray.into_iter().map(|x| x as u32).collect();

        let lcp = make_lcp_array(text, &array);
        let child = make_child_table(&lcp);

        let buckets_len = 1 << (2 * bucket_width);
        let mut prefix = vec![0; bucket_width];
        let mut buckets = Vec::with_capacity(buckets_len);

        for idx in 0..buckets_len {
            for (i, x) in prefix.iter_mut().enumerate() {
                *x = sequence::two_bit_to_code((idx >> (2 * i)) as u8);
            }
            let range = search_suffix_array(&array, &child, text, &prefix, 0, 0, text.len(), 0)
                .unwrap_or_default();
            buckets.push(range.start as u32..range.end as u32);
        }

        Self {
            array,
            child,
            buckets,
            bucket_width,
        }
    }

    pub fn search(&self, text: &[u8], query: &[u8]) -> Option<Range<usize>> {
        debug_assert!(self.bucket_width <= query.len());

        let mut idx = 0;
        for (i, x) in query[..self.bucket_width].iter().enumerate() {
            idx |= (sequence::code_to_two_bit(*x) as usize) << (2 * i);
        }

        let range = &self.buckets[idx];
        if range.start == range.end {
            return None;
        }

        let begin = range.start as usize;
        let end = range.end as usize;

        if query.len() == self.bucket_width {
            return Some(begin..end);
        }

        let store_pos = if self.child[begin] < range.end {
            begin
        } else {
            end - 1
        };

        search_suffix_array(
            &self.array,
            &self.child,
            text,
            query,
            self.bucket_width,
            begin,
            end,
            store_pos,
        )
    }

    /// Searches shortest (>= `min_len`) match of prefix of `query` to the `text`
    /// with at most `max_hits` hits.
    ///
    /// Returns `(suffix array range, match length)`.
    /// If no such match was found, `None` is returned.
    pub fn extension_search(
        &self,
        text: &[u8],
        query: &[u8],
        min_len: usize,
        max_hits: usize,
    ) -> Option<(Range<usize>, usize)> {
        debug_assert!(self.bucket_width <= min_len && min_len <= query.len());

        SEARCH_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        let mut idx = 0;
        for (i, x) in query[..self.bucket_width].iter().enumerate() {
            idx |= (sequence::code_to_two_bit(*x) as usize) << (2 * i);
        }

        let range = &self.buckets[idx];
        if range.start == range.end {
            return None;
        }

        let mut depth = self.bucket_width;
        let mut begin = range.start as usize;
        let mut end = range.end as usize;
        let mut store_pos = if self.child[begin] < range.end {
            begin
        } else {
            end - 1
        };

        while depth < min_len {
            if !do_one_step(
                &self.array,
                &self.child,
                text,
                query,
                &mut depth,
                &mut begin,
                &mut end,
                &mut store_pos,
            ) {
                return None;
            }
        }

        FOUND_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        let query_len = query.len();

        while depth < query_len && end - begin > max_hits {
            if !do_one_step(
                &self.array,
                &self.child,
                text,
                query,
                &mut depth,
                &mut begin,
                &mut end,
                &mut store_pos,
            ) {
                return None;
            }
        }

        if depth == query_len && end - begin > max_hits {
            None
        } else {
            Some((begin..end, depth))
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn search_suffix_array(
    array: &[u32],
    child: &[u32],
    text: &[u8],
    query: &[u8],
    mut depth: usize,
    mut begin: usize,
    mut end: usize,
    mut store_pos: usize,
) -> Option<Range<usize>> {
    while depth < query.len() {
        if !do_one_step(
            array,
            child,
            text,
            query,
            &mut depth,
            &mut begin,
            &mut end,
            &mut store_pos,
        ) {
            return None;
        }
    }

    Some(begin..end)
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn do_one_step(
    array: &[u32],
    child: &[u32],
    text: &[u8],
    query: &[u8],
    depth: &mut usize,
    begin: &mut usize,
    end: &mut usize,
    store_pos: &mut usize,
) -> bool {
    STEP_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

    let q = query[*depth];
    let mut b = text[array[*begin] as usize + *depth];
    if q < b {
        return false;
    }

    loop {
        let e = text[array[*end - 1] as usize + *depth];
        if q > e {
            return false;
        }

        loop {
            if b == e {
                *depth += 1;
                return true;
            }
            let mid = child[*store_pos] as usize;
            let m = text[array[mid] as usize + *depth];
            if q < m {
                *end = mid;
                *store_pos = mid - 1;
                break;
            }
            *begin = mid;
            b = m;
            *store_pos = mid;
        }
    }
}

fn make_lcp_array(text: &[u8], sa: &[u32]) -> Vec<u32> {
    let len = text.len();
    let mut inv_sa = vec![0; len];
    let mut lcp = vec![0; len];
    for i in 0..len {
        inv_sa[sa[i] as usize] = i;
    }

    // Kasai's algorithm
    let mut k = 0;
    for i in 0..len {
        if inv_sa[i] == len - 1 {
            k = 0;
            continue;
        }
        let j = sa[inv_sa[i] + 1] as usize;
        while i + k < len && j + k < len && text[i + k] == text[j + k] {
            k += 1;
        }
        lcp[inv_sa[i] + 1] = k as u32;
        if k > 0 {
            k -= 1;
        }
    }

    lcp
}

fn make_child_table(lcp: &[u32]) -> Vec<u32> {
    use std::cmp::Ordering::*;

    let mut child_table = vec![0; lcp.len()];
    let mut stack = vec![(0, lcp.len(), 0)];
    let mut min_indices = Vec::new();
    while let Some((beg, end, store_pos)) = stack.pop() {
        if end - beg < 2 {
            continue;
        }

        let mut min_lcp = u32::MAX;
        min_indices.clear();
        for (i, l) in lcp.iter().enumerate().take(end).skip(beg + 1) {
            match l.cmp(&min_lcp) {
                Less => {
                    min_lcp = *l;
                    min_indices.clear();
                    min_indices.push(i);
                }
                Equal => {
                    min_indices.push(i);
                }
                _ => {}
            }
        }

        let mid = min_indices[min_indices.len() / 2];
        child_table[store_pos] = mid as u32;

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

    #[test]
    fn lcp_array_and_child_table() {
        let text = b"gcctagccta";
        let sa = &[4, 9, 1, 6, 2, 7, 0, 5, 3, 8];
        let lcp = make_lcp_array(text, sa);
        let child = make_child_table(&lcp);
        assert_eq!(lcp, &[0, 1, 0, 4, 1, 3, 0, 5, 0, 2]);
        assert_eq!(child, &[6, 1, 4, 3, 5, 2, 8, 7, 9, 0]);
    }

    fn check_search(text: &[u8], query: &[u8], expected: Option<Vec<usize>>) {
        let text = sequence::encode(text);
        let query = sequence::encode(query);
        let sa = SuffixArray::new(&text, 1);
        let got: Option<Vec<usize>> = sa
            .search(&text, &query)
            .map(|range| range.map(|i| sa.array[i] as usize).sorted().collect());
        assert_eq!(got, expected);
    }

    #[test]
    fn search() {
        check_search(b"agctagt$", b"gct", Some(vec![1]));
        check_search(b"agctagt$", b"agctagta", None);
        check_search(b"agctagt$", b"gca", None);
        check_search(
            b"atcgggatatatggagagcttagag$",
            b"gag",
            Some(vec![13, 15, 22]),
        );
        check_search(b"agct$tgta$", b"a", Some(vec![0, 8]));
        check_search(
            b"atcgggatatatggagagcttagag$tttagagggttcttcgggatt$",
            b"gag",
            Some(vec![13, 15, 22, 30]),
        );
        check_search(
            b"atatatgca$atatatgct$atatatgga$",
            b"atat",
            Some(vec![0, 2, 10, 12, 20, 22]),
        );
        check_search(b"ca$cc$cg$", b"cc", Some(vec![3]));
        check_search(b"taat$tata$tcac$tcag$", b"taa", Some(vec![0]));
    }

    fn check_extension_search(
        text: &[u8],
        query: &[u8],
        min_len: usize,
        max_hits: usize,
        expected: Option<(Vec<usize>, usize)>,
    ) {
        let text = sequence::encode(text);
        let query = sequence::encode(query);
        let sa = SuffixArray::new(&text, 2);

        let got: Option<(Vec<_>, usize)> = sa
            .extension_search(&text, &query, min_len, max_hits)
            .map(|(range, len)| (range.map(|i| sa.array[i] as usize).sorted().collect(), len));

        assert_eq!(got, expected);
        if let Some((hits, _)) = got {
            assert!(hits.len() <= max_hits);
        }
    }

    #[test]
    fn extension_search() {
        check_extension_search(b"atatatgca$atatatgct$atatatgga$", b"atatat", 6, 2, None);
        check_extension_search(
            b"atatatgca$atatatgct$atatatgga$",
            b"atatat",
            4,
            3,
            Some((vec![0, 10, 20], 5)),
        );
        check_extension_search(b"atatatgca$atatatgct$atatatgga$", b"atatat", 4, 2, None);
        check_extension_search(
            b"atatatgca$atatatgct$atatatgga$",
            b"atatatgc",
            4,
            2,
            Some((vec![0, 10], 8)),
        );
    }
}
