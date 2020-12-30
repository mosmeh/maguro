use super::Anchor;
use itertools::Itertools;

#[derive(Debug)]
pub struct Chain {
    pub score: f64,
    pub anchors: Vec<Anchor>,
}

/// Finds optimal chains based on minimap2's formulation.
pub fn chain_anchors(
    mut anchors: Vec<Anchor>,
    max_splice_gap: usize,
    score_threshold: f64,
) -> Vec<Chain> {
    anchors.sort_by(|a, b| {
        (a.ref_pos + a.len)
            .cmp(&(b.ref_pos + b.len))
            .then_with(|| (a.query_pos + a.len).cmp(&(b.query_pos + b.len)))
    });

    let anchors: Vec<_> = coalesce_anchors(anchors.into_iter()).collect();

    let avg_seed = anchors.iter().map(|m| m.len).sum::<usize>() as f64 / anchors.len() as f64;

    let mut f = vec![0.0; anchors.len()];
    let mut ptrs = vec![0; anchors.len()];
    let mut best_score = f64::MIN;
    let mut best_chain_ends = Vec::new();

    for (i, anchor_i) in anchors.iter().enumerate() {
        let query_end_i = anchor_i.query_pos + anchor_i.len;
        let ref_end_i = anchor_i.ref_pos + anchor_i.len;

        f[i] = anchor_i.len as f64;
        ptrs[i] = i;

        for (j, anchor_j) in anchors[..i].iter().enumerate() {
            let query_end_j = anchor_j.query_pos + anchor_j.len;
            let ref_end_j = anchor_j.ref_pos + anchor_j.len;

            let query_diff = query_end_i as isize - query_end_j as isize;
            let ref_diff = ref_end_i as isize - ref_end_j as isize;

            // avoid overlap
            if query_end_j > anchor_i.query_pos || ref_end_j > anchor_i.ref_pos {
                continue;
            }

            let extension_score = f[j] + alpha(query_diff, ref_diff, anchor_i.len)
                - beta(query_diff, ref_diff, avg_seed, max_splice_gap);

            if extension_score > f[i] {
                ptrs[i] = j;
                f[i] = extension_score;
            }
        }

        if f[i] >= best_score {
            if f[i] > best_score {
                best_score = f[i];
                best_chain_ends.clear();
            }
            best_chain_ends.push(i);
        }
    }

    if best_score < score_threshold {
        return Vec::new();
    }

    // backtracking

    let mut used = vec![false; anchors.len()];
    best_chain_ends
        .into_iter()
        .filter(|end| f[*end] >= score_threshold)
        .sorted_by(|a, b| f[*b].partial_cmp(&f[*a]).unwrap())
        .filter_map(|end| {
            if used[end] {
                return None;
            }

            let mut next = ptrs[end];
            let mut curr = end;
            let mut anchor_indices = Vec::new();

            while next < curr {
                if used[curr] {
                    break;
                }
                anchor_indices.push(curr);
                used[curr] = true;
                curr = next;
                next = ptrs[next];
            }
            if !used[curr] {
                assert_eq!(curr, next);
                anchor_indices.push(curr);
            }

            let anchors = anchor_indices
                .iter()
                .rev()
                .map(|i| anchors[*i].clone())
                .collect();

            Some(Chain {
                score: f[end],
                anchors,
            })
        })
        .collect()
}

fn alpha(query_diff: isize, ref_diff: isize, len: usize) -> f64 {
    query_diff.min(ref_diff).min(len as isize) as f64
}

fn beta(query_diff: isize, ref_diff: isize, avg_seed: f64, max_splice_gap: usize) -> f64 {
    let abs_l = (query_diff - ref_diff).abs() as usize;
    if query_diff <= 0 || ref_diff <= 0 || abs_l > max_splice_gap {
        f64::INFINITY
    } else if abs_l == 0 {
        0.0
    } else {
        0.05 * avg_seed * (abs_l as f64) + 0.5 * (abs_l as f64).log2()
    }
}

fn coalesce_anchors(anchors: impl Iterator<Item = Anchor>) -> impl Iterator<Item = Anchor> {
    anchors.coalesce(|a, b| {
        let a_query_end = a.query_pos + a.len;
        let b_query_start = b.query_pos;

        let a_ref_end = a.ref_pos + a.len;
        let b_ref_start = b.ref_pos;
        let b_ref_end = b.ref_pos + b.len;

        let query_overlap = a_query_end as isize - b_query_start as isize;
        let ref_overlap = a_ref_end as isize - b_ref_start as isize;

        if ref_overlap >= 0 && (ref_overlap == query_overlap) {
            Ok(Anchor {
                ref_pos: a.ref_pos,
                query_pos: a.query_pos,
                len: a.len + b_ref_end - a_ref_end,
            })
        } else {
            Err((a, b))
        }
    })
}
