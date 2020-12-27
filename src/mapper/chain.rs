use super::Anchor;
use itertools::Itertools;

pub struct Chain {
    pub score: f64,
    pub anchors: Vec<Anchor>,
}

pub fn chain_anchors(mut anchors: Vec<Anchor>, max_splice_gap: usize) -> Vec<Chain> {
    anchors.sort_by(|a, b| (a.ref_pos + a.len).cmp(&(b.ref_pos + b.len)));

    let avg_seed = anchors.iter().map(|m| m.len).sum::<usize>() as f64 / anchors.len() as f64;

    let mut f = vec![0.0; anchors.len()];
    let mut p = vec![0; anchors.len()];
    let mut best_score = f64::MIN;
    let mut best_chain_ends = Vec::new();

    for (i, anchor_i) in anchors.iter().enumerate() {
        let query_pos_i = anchor_i.query_pos + anchor_i.len;
        let ref_pos_i = anchor_i.ref_pos + anchor_i.len;

        f[i] = anchor_i.len as f64;
        p[i] = i;

        for (j, anchor_j) in anchors[..i].iter().enumerate() {
            let query_pos_j = anchor_j.query_pos + anchor_j.len;
            let ref_pos_j = anchor_j.ref_pos + anchor_j.len;

            let query_diff = query_pos_i as isize - query_pos_j as isize;
            let ref_diff = ref_pos_i as isize - ref_pos_j as isize;
            let extension_score = f[j] + alpha(query_diff, ref_diff, anchor_i.len)
                - beta(query_diff, ref_diff, avg_seed, max_splice_gap);

            if extension_score > f[i] {
                p[i] = j;
                f[i] = extension_score;
            }
        }

        if f[i] > best_score {
            best_score = f[i];
            best_chain_ends.clear();
            best_chain_ends.push(i);
        } else if (f[i] - best_score).abs() <= f64::EPSILON {
            best_chain_ends.push(i);
        }
    }

    let mut used = vec![false; anchors.len()];
    best_chain_ends
        .iter()
        .filter_map(|end| {
            let mut ptr = p[*end];
            let mut last_ptr = *end;
            let mut anchor_indices = Vec::new();
            let mut overlap = false;

            while ptr < last_ptr {
                if used[last_ptr] {
                    overlap = true;
                }
                anchor_indices.push(last_ptr);
                used[last_ptr] = true;
                last_ptr = ptr;
                ptr = p[ptr];
            }
            if used[last_ptr] || overlap {
                return None;
            }
            anchor_indices.push(last_ptr);

            let anchors = anchor_indices
                .iter()
                .rev()
                .map(|i| anchors[*i].clone())
                .coalesce(|a, b| {
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
                .collect();

            Some(Chain {
                score: f[*end],
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
