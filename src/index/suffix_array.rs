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
    pub offsets: Vec<u32>,
    pub buckets: Vec<u32>,
    pub k: usize,
    pub f: f64,
}

impl SuffixArray {
    pub fn new(text: &[u8], k: usize, f: f64) -> Self {
        assert!(text.len() <= u32::MAX as usize + 1);

        let sa = SA::<i32>::new(text);
        let array: Vec<_> = sa.sarray.into_iter().map(|x| x as u32).collect();

        let offsets_len = 1 << (2 * k);
        let mut counts = vec![0u32; offsets_len];
        for i in 0..=(text.len() - k) {
            let seq = &text[i..i + k];
            if seq.iter().any(|x| *x == 0 || *x == sequence::DUMMY_CODE) {
                continue;
            }
            let mut idx = 0;
            for (j, x) in seq.iter().enumerate() {
                idx |= (sequence::code_to_two_bit(*x) as usize) << (2 * j);
            }
            counts[idx] += 1;
        }

        let mut buckets_len = 0;
        let mut offsets = Vec::new();
        for i in 0..offsets_len {
            let w = ((counts[i] as f64 * f).log(4.0).max(0.0) as usize).min(31);
            offsets.push(buckets_len as u32);
            buckets_len += 1 << (2 * w);
        }
        offsets.push(buckets_len as u32);

        let mut ssa = Vec::new();
        let mut buckets = vec![u32::MAX; buckets_len];
        let mut prev_bucket = 0;
        for (i, s) in array.into_iter().enumerate() {
            if s as usize + k > text.len() {
                continue;
            }
            let seq = &text[s as usize..s as usize + k];
            if seq.iter().any(|x| *x == 0 || *x == sequence::DUMMY_CODE) {
                continue;
            }
            let mut idx = 0;
            for (j, x) in seq.iter().rev().enumerate() {
                idx |= (sequence::code_to_two_bit(*x) as usize) << (2 * j);
            }
            let w = ((offsets[idx + 1] - offsets[idx]).trailing_zeros() / 2) as usize;
            if s as usize + w > text.len() {
                continue;
            }
            let seq2 = &text[s as usize + k..][..w];
            if seq2.iter().any(|x| *x == 0 || *x == sequence::DUMMY_CODE) {
                continue;
            }
            let mut idx2 = 0;
            for (j, x) in seq2.iter().rev().enumerate() {
                idx2 |= (sequence::code_to_two_bit(*x) as usize) << (2 * j);
            }
            assert!(idx2 < (1 << (2 * w)));
            let j = offsets[idx] as usize + idx2;
            /*println!(
                "{} {} {}",
                idx,
                String::from_utf8(crate::sequence::decode(
                    &text[s as usize..s as usize + k + ws[idx]]
                ))
                .unwrap(),
                j
            );*/
            assert!(prev_bucket <= j, "{} {}", prev_bucket, j);
            if buckets[j] == u32::MAX {
                buckets[j] = ssa.len() as u32;
            }
            prev_bucket = j;
            ssa.push(s);
        }
        buckets.push(ssa.len() as u32);

        for i in (0..buckets_len).rev() {
            if buckets[i] == u32::MAX {
                buckets[i] = buckets[i + 1];
            }
        }

        for i in 0..offsets_len {
            assert!(offsets[i] <= offsets[i + 1]);
        }
        for i in 0..buckets_len {
            assert!(buckets[i] <= buckets[i + 1]);
        }

        Self {
            array: ssa,
            offsets,
            buckets,
            k,
            f,
        }
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
        debug_assert!(self.k <= min_len && min_len <= query.len());

        SEARCH_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        let mut idx = 0;
        for (i, x) in query[..self.k].iter().rev().enumerate() {
            idx |= (sequence::code_to_two_bit(*x) as usize) << (2 * i);
        }

        let bucket_begin = self.offsets[idx] as usize;
        let bucket_end = self.offsets[idx + 1] as usize;
        if bucket_begin == bucket_end {
            return None;
        }

        let w = ((bucket_end - bucket_begin).trailing_zeros() / 2) as usize;
        let mut idx2 = 0;
        for (i, x) in query[self.k..self.k + w].iter().rev().enumerate() {
            idx2 |= (sequence::code_to_two_bit(*x) as usize) << (2 * i);
        }

        let mut begin = self.buckets[bucket_begin + idx2] as usize;
        let mut end = self.buckets[bucket_begin + idx2 + 1] as usize;
        if begin == end {
            return None;
        }

        unsafe {
            equal_range(
                &self.array,
                text.as_ptr().add(self.k + w),
                query.as_ptr().add(self.k + w),
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
                    &self.array,
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
