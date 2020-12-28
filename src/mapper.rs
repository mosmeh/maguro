pub mod align;
mod chain;

use crate::index::{Index, SequenceId};
use align::{Aligner, AlignmentConfig};
use chain::Chain;
use itertools::Itertools;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Anchor {
    query_pos: usize,
    ref_pos: usize,
    len: usize,
}

pub struct Mapping {
    pub seq_id: SequenceId,
    pub pos: usize,
}

pub struct Mapper<'a> {
    index: &'a Index,
    k: usize,
    consensus_fraction: f64,
    coverage_score_ratio: f64,
    max_splice_gap: usize,
    min_score_fraction: f64,
    aligner: Aligner,
}

impl Mapper<'_> {
    pub fn map(&self, query: &[u8]) -> Vec<Mapping> {
        if query.len() < self.k {
            unimplemented!()
        }

        let mut ref_anchors: HashMap<SequenceId, Vec<Anchor>> = HashMap::new();

        for kmer_pos in 0..=(query.len() - self.k) {
            let range = self
                .index
                .sa
                .search(&self.index.seq, &query[kmer_pos..kmer_pos + self.k]);
            if let Some(range) = range {
                for i in range {
                    let pos = self.index.sa.array[i];
                    let id = self.index.seq_id_from_pos(pos);

                    ref_anchors
                        .entry(id)
                        .and_modify(|anchors| {
                            anchors.push(Anchor {
                                query_pos: kmer_pos,
                                ref_pos: pos,
                                len: self.k,
                            })
                        })
                        .or_insert_with(|| {
                            vec![Anchor {
                                query_pos: kmer_pos,
                                ref_pos: pos,
                                len: self.k,
                            }]
                        });
                }
            }
        }

        let mut ref_chains: HashMap<SequenceId, Vec<Chain>> = HashMap::new();
        let mut global_best_score = f64::MIN;

        for (seq_id, anchors) in ref_anchors.into_iter() {
            let min_accepted_score_global = global_best_score * self.consensus_fraction;
            let mut chains =
                chain::chain_anchors(anchors, self.max_splice_gap, min_accepted_score_global);
            if chains.is_empty() {
                continue;
            }

            let best_score = chains
                .iter()
                .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
                .unwrap()
                .score;
            global_best_score = global_best_score.max(best_score);

            let min_accepted_score_per_ref = best_score * self.coverage_score_ratio;
            chains.retain(|chain| chain.score >= min_accepted_score_per_ref);

            if !chains.is_empty() {
                ref_chains.insert(seq_id, chains);
            }
        }

        let min_accepted_chain_score = global_best_score * self.consensus_fraction;

        let match_score = self.aligner.config().match_score as i32;
        let max_align_score = query.len() as i32 * match_score;
        let min_accepted_align_score = (max_align_score as f64 * self.min_score_fraction) as i32;

        let mut rev_query = query.to_owned();
        rev_query.reverse();

        let mut mappings = Vec::new();
        for (seq_id, chains) in ref_chains.iter() {
            let seq_range = self.index.seq_range(*seq_id);
            let chains = chains
                .iter()
                .filter(|chain| chain.score >= min_accepted_chain_score);

            let mut rev_seq = self.index.seq[seq_range.start..seq_range.end].to_owned();
            rev_seq.reverse();

            for chain in chains {
                let anchors = &chain.anchors;
                let mut align_score = 0;

                // exact matches (anchors)
                align_score +=
                    match_score * anchors.iter().map(|anchor| anchor.len as i32).sum::<i32>();

                // between anchors
                align_score += anchors
                    .iter()
                    .tuple_windows()
                    .map(|(a, b)| {
                        self.aligner.global_align(
                            &query[a.query_pos + a.len..b.query_pos],
                            &self.index.seq[a.ref_pos + a.len..b.ref_pos],
                        )
                    })
                    .sum::<i32>();

                let first_anchor = &anchors[0];
                let last_anchor = &anchors[anchors.len() - 1];
                let query_anchor_range =
                    first_anchor.query_pos..(last_anchor.query_pos + last_anchor.len);
                let ref_anchor_range =
                    first_anchor.ref_pos..(last_anchor.ref_pos + last_anchor.len);

                // left extension
                if query_anchor_range.start > 0 {
                    let target = if ref_anchor_range.start > seq_range.start {
                        &rev_seq[rev_seq.len() - (ref_anchor_range.start - seq_range.start) + 1..]
                    } else {
                        &[]
                    };
                    align_score += self.aligner.extension_align(
                        &rev_query[rev_query.len() - query_anchor_range.start + 1..],
                        target,
                    );
                }

                // right extension
                if query_anchor_range.end < query.len() {
                    align_score += self.aligner.extension_align(
                        &query[query_anchor_range.end..],
                        &self.index.seq[ref_anchor_range.end..seq_range.end],
                    );
                }

                if align_score >= min_accepted_align_score {
                    mappings.push(Mapping {
                        seq_id: *seq_id,
                        pos: (ref_anchor_range.start - seq_range.start)
                            .saturating_sub(query_anchor_range.start),
                    })
                }
            }
        }

        mappings
    }
}

pub struct MapperBuilder<'a> {
    index: &'a Index,
    k: usize,
    consensus_fraction: f64,
    coverage_score_ratio: f64,
    max_splice_gap: usize,
    min_score_fraction: f64,
    alignment_config: AlignmentConfig,
}

impl<'a> MapperBuilder<'a> {
    pub fn new(index: &'a Index) -> Self {
        Self {
            index,
            k: 15,
            consensus_fraction: 0.65,
            coverage_score_ratio: 0.6,
            max_splice_gap: 100,
            min_score_fraction: 0.65,
            alignment_config: Default::default(),
        }
    }

    pub fn k(&mut self, k: usize) -> &mut Self {
        self.k = k;
        self
    }

    pub fn consensus_fraction(&mut self, consensus_fraction: f64) -> &mut Self {
        self.consensus_fraction = consensus_fraction;
        self
    }

    pub fn coverage_score_ratio(&mut self, coverage_score_ratio: f64) -> &mut Self {
        self.coverage_score_ratio = coverage_score_ratio;
        self
    }

    pub fn max_splice_gap(&mut self, max_splice_gap: usize) -> &mut Self {
        self.max_splice_gap = max_splice_gap;
        self
    }

    pub fn min_score_fraction(&mut self, min_score_fraction: f64) -> &mut Self {
        self.min_score_fraction = min_score_fraction;
        self
    }

    pub fn alignment_config(&mut self, alignment_config: AlignmentConfig) -> &mut Self {
        self.alignment_config = alignment_config;
        self
    }

    pub fn build(&self) -> Mapper<'a> {
        Mapper {
            index: self.index,
            k: self.k,
            consensus_fraction: self.consensus_fraction,
            coverage_score_ratio: self.coverage_score_ratio,
            max_splice_gap: self.max_splice_gap,
            min_score_fraction: self.min_score_fraction,
            aligner: Aligner::new(self.alignment_config.clone()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{alphabet, index::IndexBuilder};
    use itertools::Itertools;

    fn check_map(seqs: &[(&[u8], &[u8])], query: &[u8], expected: &[&[u8]]) {
        let mut fasta: Vec<u8> = Vec::new();
        for (name, seq) in seqs {
            fasta.push(b'>');
            fasta.extend(*name);
            fasta.push(b'\n');
            fasta.extend(*seq);
            fasta.push(b'\n');
        }

        let cursor = std::io::Cursor::new(&fasta);
        let index = IndexBuilder::new(cursor).build();
        let mapper = MapperBuilder::new(&index)
            .k(3)
            .min_score_fraction(0.3)
            .build();

        let query = alphabet::encode_seq(&query);

        let got: Vec<_> = mapper
            .map(&query)
            .iter()
            .map(|mapping| mapping.seq_id)
            .unique()
            .map(|seq_id| index.seq_name(seq_id))
            .sorted()
            .collect();
        let expected: Vec<_> = expected.iter().sorted().copied().collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn map() {
        check_map(&[(b"foo", b"agctagt")], b"gct", &[b"foo"]);
        check_map(&[(b"foo", b"agctagt")], b"gctagt", &[b"foo"]);
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"ggatcgatggagctctt",
            &[b"foo"],
        );
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"tcgatataggagcgctta",
            &[b"foo"],
        );
        check_map(&[(b"foo", b"agct"), (b"bar", b"tgta")], b"gct", &[b"foo"]);
        check_map(
            &[
                (b"foo", b"atcgggatatatggagagcttagag"),
                (b"bar", b"tatggagggagccttag"),
            ],
            b"atggagagcttag",
            &[b"foo", b"bar"],
        );
    }
}
