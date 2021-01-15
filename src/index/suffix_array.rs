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
    pub prefix_len: usize,
    pub table1: Vec<u32>,
    pub table2: Vec<Range<u32>>,
}

impl SuffixArray {
    pub fn new(text: &[u8], prefix_len: usize, f: f64) -> Self {
        assert!(text.len() <= u32::MAX as usize + 1);
        assert!(prefix_len < 32);

        let sa = SA::<i32>::new(text);
        let array: Vec<_> = sa.sarray.into_iter().map(|x| x as u32).collect();

        let lcp = make_lcp_array(text, &array);
        let child = make_child_table(&lcp);
        eprintln!("Built SA");

        let table1_len = 1 << (2 * prefix_len);
        let mut prefix = vec![0; prefix_len];
        let mut table1 = Vec::with_capacity(table1_len);
        let mut table2 = Vec::new();

        let mut num_found = 0;
        let mut hist_w = std::collections::HashMap::new();
        let mut min_occ = usize::MAX;
        let mut max_occ = usize::MIN;
        let mut sum_occ = 0;
        let mut theo_reduction_vl = 0.0;
        let mut real_reduction_vl = 0.0;
        let mut non_empty_entries = 0;
        for idx in 0..table1_len {
            for (i, x) in prefix.iter_mut().enumerate() {
                *x = sequence::two_bit_to_code((idx >> (2 * i)) as u8);
            }
            if let Some(range) =
                search_suffix_array(&array, &child, text, &prefix, 0, 0, text.len(), 0)
            {
                num_found += 1;
                let occ = range.end - range.start;
                min_occ = min_occ.min(occ);
                max_occ = max_occ.max(occ);
                sum_occ += occ;
                theo_reduction_vl += occ as f64 * (prefix_len as f64 + (occ as f64).log(4.0));

                let freq = (range.end - range.start) as f64 * f;
                let w = (freq.log(4.0).max(0.0) as usize).min(31);
                real_reduction_vl += occ as f64 * (prefix_len as f64 + w as f64);
                hist_w.entry(w).and_modify(|x| *x += 1).or_insert(1);
                table1.push(table2.len() as u32);
                if w == 0 {
                    table2.push(range.start as u32..range.end as u32);
                    non_empty_entries += 1;
                } else {
                    let mut body = vec![0; prefix_len + w];
                    for offset in 0..(1 << (2 * w)) {
                        for (i, x) in body.iter_mut().skip(prefix_len).enumerate() {
                            *x = sequence::two_bit_to_code((offset >> (2 * i)) as u8);
                        }
                        let range2 = search_suffix_array(
                            &array,
                            &child,
                            text,
                            &body,
                            prefix_len,
                            range.start,
                            range.end,
                            if child[range.start as usize] < range.end as u32 {
                                range.start
                            } else {
                                range.end - 1
                            },
                        )
                        .unwrap_or_default();
                        table2.push(range2.start as u32..range2.end as u32);
                        if range2.end > range2.start {
                            non_empty_entries += 1;
                        }
                    }
                }
            } else {
                min_occ = 0;
                hist_w.entry(0).and_modify(|x| *x += 1).or_insert(1);
                table1.push(table2.len() as u32);
            }
        }

        eprintln!("k={}, f={}", prefix_len, f);
        eprintln!(
            "Prefix counts: min={} max={} avg={:.2}",
            min_occ,
            max_occ,
            sum_occ as f64 / table1_len as f64
        );
        eprintln!("Reduction:");
        theo_reduction_vl /= sum_occ as f64;
        eprintln!(
            "Theo ord={:.4} vl={:.4}",
            (text.len() as f64).log(4.0),
            theo_reduction_vl
        );
        real_reduction_vl /= sum_occ as f64;
        eprintln!(
            "Real ord={:.4} vl={:.4}",
            (text.len() as f64 * f).log(4.0).floor(),
            real_reduction_vl
        );
        eprintln!(
            "Found {} / {} prefixes ({:.2}%)",
            num_found,
            table1.len(),
            num_found as f64 * 100.0 / table1.len() as f64
        );
        eprintln!("Table 2 summary:");
        use itertools::Itertools;
        for (w, count) in hist_w.iter().sorted_by(|(a, _), (b, _)| a.cmp(&b)) {
            eprintln!("w={}\tcount={}", w, count);
        }
        eprintln!(
            "Non-empty: {} / {} ({:.2}%)",
            non_empty_entries,
            table2.len(),
            non_empty_entries as f64 * 100.0 / table2.len() as f64
        );
        eprintln!(
            "Table1 len = {}\t{}B\t{:.2}KiB",
            table1.len(),
            table1.len() * 4,
            table1.len() as f64 * 4.0 / 1024.0
        );
        eprintln!(
            "Table2 len = {}\t{}B\t{:.2}MiB",
            table2.len(),
            table2.len() * 8,
            table2.len() as f64 * 8.0 / 1024.0 / 1024.0
        );

        table1.push(table2.len() as u32);

        Self {
            array,
            child,
            prefix_len,
            table1,
            table2,
        }
    }

    pub fn search(&self, text: &[u8], query: &[u8]) -> Option<Range<usize>> {
        debug_assert!(self.prefix_len <= query.len());

        let mut start_idx = 0;
        for (i, x) in query[..self.prefix_len].iter().enumerate() {
            start_idx |= (sequence::code_to_two_bit(*x) as usize) << (2 * i);
        }

        let a = &self.table1[start_idx];
        let b = &self.table1[start_idx + 1];
        if a == b {
            return None;
        }

        let w = ((b - a).trailing_zeros() / 2) as usize;
        debug_assert!(self.prefix_len + w <= query.len());

        let mut offset = 0;
        for (i, x) in query[self.prefix_len..][..w].iter().enumerate() {
            offset |= (sequence::code_to_two_bit(*x) as usize) << (2 * i);
        }

        let range = &self.table2[*a as usize..][offset];
        if range.start == range.end {
            return None;
        }

        let depth = self.prefix_len + w;
        let begin = range.start as usize;
        let end = range.end as usize;

        if query.len() == depth {
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
            depth,
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
        debug_assert!(self.prefix_len <= min_len && min_len <= query.len());

        SEARCH_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        let mut idx = 0;
        for (i, x) in query[..self.prefix_len].iter().enumerate() {
            idx |= (sequence::code_to_two_bit(*x) as usize) << (2 * i);
        }

        let a = &self.table1[idx];
        let b = &self.table1[idx + 1];
        if a == b {
            return None;
        }

        let w = ((b - a).trailing_zeros() / 2) as usize;
        debug_assert!(self.prefix_len + w <= min_len);

        let mut offset = 0;
        for (i, x) in query[self.prefix_len..][..w].iter().enumerate() {
            offset |= (sequence::code_to_two_bit(*x) as usize) << (2 * i);
        }

        let range = &self.table2[*a as usize..][offset];
        if range.start == range.end {
            return None;
        }

        let mut depth = self.prefix_len + w;
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
        let sa = SuffixArray::new(&text, 1, 2.0);
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
        let sa = SuffixArray::new(&text, 2, 2.0);

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
