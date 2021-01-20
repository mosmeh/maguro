use serde::{Deserialize, Serialize};
use std::ops::Range;
use sufsort_rs::sufsort::SA;
use xxhash_rust::xxh32::xxh32;

pub static STEP_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub static FOUND_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub static SEARCH_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

#[derive(Serialize, Deserialize)]
pub struct SuffixArray {
    pub ssa: Vec<u32>,
    pub child: Vec<u32>,
    pub offsets: Vec<u32>,
    pub skip: Vec<u32>,
    pub k: usize,
    pub mask: usize,
}

impl SuffixArray {
    pub fn new(text: &[u8], k: usize, bits: usize) -> Self {
        assert!(text.len() <= u32::MAX as usize + 1);

        let sa = SA::<i32>::new(text);
        let array: Vec<_> = sa.sarray.into_iter().map(|x| x as u32).collect();
        let lcp = make_lcp_array(text, &array);

        let hashtable_len = 1 << bits;
        let mask = hashtable_len - 1;
        let mut counts = vec![0u32; hashtable_len];
        for i in 0..=(text.len() - k) {
            let seq = &text[i..i + k];
            if seq
                .iter()
                .any(|x| *x == 0 || *x == crate::sequence::DUMMY_CODE)
            {
                continue;
            }
            counts[xxh32(&seq, 0) as usize & mask] += 1;
        }

        let mut cum_sum = 0;
        for count in counts.iter_mut() {
            let x = *count;
            *count = cum_sum;
            cum_sum += x;
        }
        counts.push(cum_sum);

        let offsets = counts;

        let mut pos = offsets.clone();
        let mut ssa = vec![0; array.len()];
        let mut slcp = vec![0; array.len()];
        let mut skip = vec![k as u32; hashtable_len];
        let mut idx = usize::MAX;
        let len = text.len() as u32;
        for (i, s) in array.iter().enumerate() {
            if *s as usize + k > text.len() {
                continue;
            }
            let seq = &text[*s as usize..][..k];
            if seq
                .iter()
                .any(|x| *x == 0 || *x == crate::sequence::DUMMY_CODE)
            {
                continue;
            }
            if lcp[i] >= k as u32 {
                if idx == usize::MAX {
                    idx = xxh32(&seq, 0) as usize & mask;
                }
                let p = pos[idx] as usize;
                ssa[p] = *s;
                slcp[p] = lcp[i];
                if skip[idx] > lcp[i] {
                    skip[idx] = lcp[i];
                }
                pos[idx] += 1;
            } else {
                idx = xxh32(&seq, 0) as usize & mask;
                let p = pos[idx] as usize;
                ssa[p] = *s;
                slcp[p] = if p as u32 > offsets[idx] {
                    let prev = ssa[p - 1];
                    let mut l = 0;
                    while *s + l < len
                        && prev + l < len
                        && text[(*s + l) as usize] == text[(prev + l) as usize]
                    {
                        l += 1;
                    }
                    if skip[idx] > l {
                        skip[idx] = l;
                    }
                    l
                } else {
                    0
                };
                pos[idx] += 1;
            }
        }

        let mut child = vec![0; ssa.len()];
        for i in 0..(offsets.len() - 1) {
            let begin = offsets[i] as usize;
            let end = offsets[i + 1] as usize;
            if begin < end {
                make_child_table(&slcp, &mut child, begin, end);
            }
        }

        Self {
            ssa,
            child,
            offsets,
            skip,
            k,
            mask,
        }
    }

    pub fn search(&self, text: &[u8], query: &[u8]) -> Option<Range<usize>> {
        let hash = xxh32(&query[..self.k], 0) as usize & self.mask;

        let begin = self.offsets[hash] as usize;
        let end = self.offsets[hash + 1] as usize;

        if begin == end {
            return None;
        }

        search_suffix_array(&self.ssa, &self.child, text, query, 0, begin, end, begin)
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
        SEARCH_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        let hash = xxh32(&query[..self.k], 0) as usize & self.mask;
        let mut begin = self.offsets[hash] as usize;
        let mut end = self.offsets[hash + 1] as usize;
        if begin == end {
            return None;
        }

        unsafe {
            equal_range(
                &self.ssa,
                text.as_ptr(),
                query.as_ptr(),
                query.as_ptr().add(min_len),
                &mut begin,
                &mut end,
            );
        }
        if begin == end {
            return None;
        }

        FOUND_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        let mut depth = min_len;
        let query_len = query.len();
        while depth < query_len && end - begin > max_hits {
            unsafe {
                equal_range(
                    &self.ssa,
                    text.as_ptr().add(depth),
                    query.as_ptr().add(depth),
                    query.as_ptr().add(depth + 1),
                    &mut begin,
                    &mut end,
                );
            }

            if begin == end {
                return None;
            }
            depth += 1;
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

unsafe fn equal_range(
    sa: &[u32],
    text_base: *const u8,
    query_begin: *const u8,
    query_end: *const u8,
    begin: &mut usize,
    end: &mut usize,
) {
    let mut q_begin = query_begin;
    let mut q_end = q_begin;
    let mut t_begin = text_base;
    let mut t_end = t_begin;

    while begin < end {
        let mid = *begin + (*end as isize - *begin as isize) as usize / 2;
        let offset = sa[mid] as usize;
        let mut q;
        let mut t;
        if q_begin < q_end {
            q = q_begin;
            t = t_begin.add(offset);
        } else {
            q = q_end;
            t = t_end.add(offset);
        }
        let mut x;
        let mut y;
        loop {
            x = *t;
            y = *q;
            if x != y {
                break;
            };
            q = q.add(1);
            if q == query_end {
                *begin = lower_bound(sa, t_begin, q_begin, query_end, *begin, mid);
                *end = upper_bound(sa, t_end, q_end, query_end, mid + 1, *end);
                return;
            }
            t = t.add(1);
        }
        if x < y {
            *begin = mid + 1;
            q_begin = q;
            t_begin = t.sub(offset);
        } else {
            *end = mid;
            q_end = q;
            t_end = t.sub(offset);
        }
    }
}

unsafe fn lower_bound(
    sa: &[u32],
    mut text_base: *const u8,
    mut query_begin: *const u8,
    query_end: *const u8,
    mut begin: usize,
    mut end: usize,
) -> usize {
    while begin < end {
        let mid = begin + (end - begin) / 2;
        let offset = sa[mid];
        let mut t = text_base.offset(offset as isize);
        let mut q = query_begin;
        loop {
            if *t < *q {
                begin = mid + 1;
                query_begin = q;
                text_base = t.sub(offset as usize);
                break;
            }
            q = q.add(1);
            if q == query_end {
                end = mid;
                break;
            }
            t = t.add(1);
        }
    }
    begin
}

unsafe fn upper_bound(
    sa: &[u32],
    mut text_base: *const u8,
    mut query_begin: *const u8,
    query_end: *const u8,
    mut begin: usize,
    mut end: usize,
) -> usize {
    while begin < end {
        let mid = begin + (end - begin) / 2;
        let offset = sa[mid];
        let mut t = text_base.offset(offset as isize);
        let mut q = query_begin;
        loop {
            if *t > *q {
                end = mid;
                query_begin = q;
                text_base = t.sub(offset as usize);
                break;
            }
            q = q.add(1);
            if q == query_end {
                begin = mid + 1;
                break;
            }
            t = t.add(1);
        }
    }
    end
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

fn make_child_table(lcp: &[u32], child_table: &mut [u32], begin: usize, end: usize) {
    use std::cmp::Ordering::*;

    let mut stack = vec![(begin, end, begin)];
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
        let mut child = vec![0; sa.len()];
        make_child_table(&lcp, &mut child, 0, sa.len());
        assert_eq!(lcp, &[0, 1, 0, 4, 1, 3, 0, 5, 0, 2]);
        assert_eq!(child, &[6, 1, 4, 3, 5, 2, 8, 7, 9, 0]);
    }

    fn check_search(text: &[u8], query: &[u8], expected: Option<Vec<usize>>) {
        let text = sequence::encode(text);
        let query = sequence::encode(query);
        let sa = SuffixArray::new(&text, 1, 8);
        let got: Option<Vec<usize>> = sa
            .search(&text, &query)
            .map(|range| range.map(|i| sa.ssa[i] as usize).sorted().collect());
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
        let sa = SuffixArray::new(&text, 2, 8);

        let got: Option<(Vec<_>, usize)> = sa
            .extension_search(&text, &query, min_len, max_hits)
            .map(|(range, len)| (range.map(|i| sa.ssa[i] as usize).sorted().collect(), len));

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
