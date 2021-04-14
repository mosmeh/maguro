use crate::sequence;
use serde::{Deserialize, Serialize};
use std::ops::Range;
use sufsort_rs::sufsort::SA;

pub static STEP_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub static FOUND_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub static SEARCH_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

fn hash_64(mut key: u64, mask: u64) -> u64 {
    key = (!key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

#[derive(Serialize, Deserialize)]
pub struct SuffixArray {
    pub ssa: Vec<u32>,
    pub offsets: Vec<u32>,
    pub k: usize,
    pub l: usize,
}

impl SuffixArray {
    pub fn new(text: &[u8], k: usize, l: usize) -> Self {
        assert!(text.len() <= u32::MAX as usize + 1);
        assert!(k < 31);
        assert!(l < 16);
        assert!(l < k);

        let sa = SA::<i32>::new(text);
        let array: Vec<_> = sa.sarray.into_iter().map(|x| x as u32).collect();

        let offsets_len = 1 << (2 * l);
        let mut left_to_indices = vec![sorted_list::SortedList::new(); offsets_len];
        for s in array {
            if s as usize + k > text.len() {
                continue;
            }
            let seq = &text[s as usize..][..k];
            if seq
                .iter()
                .any(|x| *x == 0 || *x == crate::sequence::DUMMY_CODE)
            {
                continue;
            }

            let mut bases = 0;
            for (j, x) in seq[..k].iter().enumerate() {
                bases |= (sequence::code_to_two_bit(*x) as u64) << (2 * j);
            }
            //assert!(bases < (1 << (2 * k)));
            bases = hash_64(bases, (1 << (2 * k)) - 1);

            let left = (bases & ((1 << (2 * l)) - 1)) as u32;
            let right = (bases >> (2 * l)) as u32;
            //assert!(right as u64 == right as u64 & ((1 << (2 * (k - l))) - 1));
            left_to_indices[left as usize].insert(right, s);
        }

        let mut left_to_right_counts = Vec::new();
        for indices in &left_to_indices {
            let mut prev = None;
            let mut count = 0;
            for (i, _) in indices.iter() {
                if prev.is_none() || prev.unwrap() != i {
                    count += 1;
                    prev = Some(i);
                }
            }
            left_to_right_counts.push(count);
        }

        let mut offsets = Vec::new();
        let mut offset_start = 0;
        for i in 0..offsets_len {
            offsets.push(offset_start);
            offset_start += left_to_indices[i].len() as u32 + left_to_right_counts[i] * 2;
        }
        offsets.push(offset_start);

        let mut ssa = vec![0u32; offset_start as usize];
        for i in 0..offsets_len {
            let right_count = left_to_right_counts[i] as usize;
            if right_count == 0 {
                continue;
            }

            let mut prev = *left_to_indices[i].iter().next().unwrap().0;
            let mut z = offsets[i] as usize;
            let mut pos = z + 2 * right_count;

            ssa[z] = pos as u32;
            ssa[z + right_count] = prev;
            z += 1;

            for (right, s) in left_to_indices[i].iter() {
                if prev != *right {
                    ssa[z] = pos as u32;
                    ssa[z + right_count] = *right;
                    z += 1;
                    prev = *right;
                }
                ssa[pos] = *s;
                pos += 1;
            }

            /*assert!(pos == offsets[i + 1] as usize);
            assert!((ssa[offsets[i] as usize] - offsets[i]) % 2 == 0);*/
        }

        /*for i in 0..offsets_len {
            if offsets[i] == offsets[i + 1] {
                continue;
            }
            let i_begin = ssa[offsets[i] as usize] as usize;
            let i_end = offsets[i + 1] as usize;
            for j in &ssa[i_begin..i_end] {
                assert!((*j as usize) < text.len());
            }
            let j_begin = offsets[i] as usize;
            assert!((i_begin - j_begin) % 2 == 0);
            let j_end = (i_begin as usize - j_begin) / 2 + j_begin;
            for j in &ssa[j_begin..j_end] {
                assert!((*j as usize) < i_end);
            }
        }*/

        Self { ssa, offsets, k, l }
    }

    /*pub fn search(&self, text: &[u8], query: &[u8]) -> Option<Range<usize>> {
        let hash = xxh32(&query[..self.k], 0) as usize & self.mask;

        let begin = self.offsets[hash] as usize;
        let end = self.offsets[hash + 1] as usize;

        if begin == end {
            return None;
        }

        search_suffix_array(&self.ssa, &self.child, text, query, 0, begin, end, begin)
    }*/

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
        //assert!(query.len() >= self.k && min_len >= self.k);
        SEARCH_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        let mut bases = 0;
        for (j, x) in query[..self.k].iter().enumerate() {
            bases |= (sequence::code_to_two_bit(*x) as u64) << (2 * j);
        }
        bases = hash_64(bases, (1 << (2 * self.k)) - 1);

        let left = (bases & ((1 << (2 * self.l)) - 1)) as u32;
        let section_begin = self.offsets[left as usize];
        let section_end = self.offsets[left as usize + 1];
        if section_begin == section_end {
            return None;
        }

        let head_begin = section_begin;
        let head_end = self.ssa[section_begin as usize];
        let num_rights = (head_end - head_begin) / 2;
        let right_begin = (head_begin + num_rights) as usize;
        let right_end = head_end as usize;

        let right = (bases >> (2 * self.l)) as u32;
        let idx = if let Ok(i) = self.ssa[right_begin..right_end].binary_search(&right) {
            i
        } else {
            return None;
        };

        let mut begin = self.ssa[head_begin as usize + idx] as usize;
        let mut end = if head_begin as usize + idx + 1 == right_begin {
            section_end
        } else {
            self.ssa[head_begin as usize + idx + 1]
        } as usize;

        /*println!(
            "{}",
            String::from_utf8(crate::sequence::decode(&query)).unwrap()
        );
        println!("{} {}", begin, end);
        for i in &self.ssa[begin..end] {
            println!(
                "{}",
                String::from_utf8(crate::sequence::decode(
                    &text[*i as usize..*i as usize + 30]
                ))
                .unwrap()
            );
        }*/

        unsafe {
            equal_range(
                &self.ssa,
                text.as_ptr().add(self.k),
                query.as_ptr().add(self.k),
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
