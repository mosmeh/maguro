pub mod align;
mod chain;

use crate::{
    index::{Index, SequenceId},
    sequence,
};
use align::{Aligner, AlignmentConfig};
use chain::Chain;
use itertools::Itertools;
use std::collections::HashMap;

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub fn is_forward(self) -> bool {
        matches!(self, Strand::Forward)
    }

    pub fn is_reverse(self) -> bool {
        matches!(self, Strand::Reverse)
    }
}

#[derive(Clone, Debug)]
pub struct Anchor {
    query_pos: usize,
    ref_pos: usize,
    len: usize,
}

pub struct Mapping {
    pub seq_id: SequenceId,
    pub pos: usize,
    pub strand: Strand,
    pub score: i32,
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
    pub fn map(&mut self, query: &[u8]) -> Vec<Mapping> {
        assert!(query.len() >= self.k);

        let rc_query = sequence::reverse_complement(&query);

        // seeding

        let ref_to_anchors = self.search_anchors(query, &rc_query);
        if ref_to_anchors.is_empty() {
            return Vec::new();
        }

        // chaining

        let mut ref_to_chains: HashMap<(SequenceId, Strand), Vec<Chain>> = HashMap::new();
        let mut best_chain_score = f64::MIN;
        let mut chain_score_threshold = f64::MIN;

        for ((seq_id, strand), anchors) in ref_to_anchors.into_iter() {
            let mut chains =
                chain::chain_anchors(anchors, self.max_splice_gap, chain_score_threshold);
            if chains.is_empty() {
                continue;
            }

            let ref_best_score = chains
                .iter()
                .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
                .unwrap()
                .score;
            best_chain_score = best_chain_score.max(ref_best_score);
            chain_score_threshold = best_chain_score * self.consensus_fraction;

            let ref_chain_score_threshold =
                chain_score_threshold.max(ref_best_score * self.coverage_score_ratio);
            chains.retain(|chain| chain.score >= ref_chain_score_threshold);

            if !chains.is_empty() {
                ref_to_chains.insert((seq_id, strand), chains);
            }
        }

        if ref_to_chains.is_empty() {
            return Vec::new();
        }

        // base-to-base alignment

        let rev_query = sequence::reverse(query);
        let compl_query = sequence::complement(query);

        let match_score = self.aligner.config().match_score as i32;
        let max_align_score = query.len() as i32 * match_score;
        let align_score_threshold = (max_align_score as f64 * self.min_score_fraction) as i32;

        let mut mappings = Vec::new();

        for ((seq_id, strand), chains) in ref_to_chains.iter() {
            let (align_query, rev_align_query): (&[u8], &[u8]) = match strand {
                Strand::Forward => (query, &rev_query),
                Strand::Reverse => (&rc_query, &compl_query),
            };

            let chains = chains
                .iter()
                .filter(|chain| chain.score >= chain_score_threshold);

            if let Some(mapping) = self.calc_best_mapping(
                *seq_id,
                *strand,
                align_query,
                rev_align_query,
                chains,
                align_score_threshold,
            ) {
                if mapping.score >= align_score_threshold {
                    mappings.push(mapping);
                }
            }
        }

        mappings
    }

    fn search_anchors(
        &self,
        query: &[u8],
        rc_query: &[u8],
    ) -> HashMap<(SequenceId, Strand), Vec<Anchor>> {
        let mut ref_to_anchors: HashMap<(SequenceId, Strand), Vec<Anchor>> = HashMap::new();

        let mut seed = |query: &[u8], strand| {
            for seed_pos in 0..=(query.len() - self.k) {
                let range = self
                    .index
                    .sa
                    .search(&self.index.seq, &query[seed_pos..seed_pos + self.k]);
                if let Some(range) = range {
                    for i in range {
                        let pos = self.index.sa.array[i];
                        let id = self.index.seq_id_from_pos(pos);

                        ref_to_anchors
                            .entry((id, strand))
                            .and_modify(|anchors| {
                                anchors.push(Anchor {
                                    query_pos: seed_pos,
                                    ref_pos: pos,
                                    len: self.k,
                                })
                            })
                            .or_insert_with(|| {
                                vec![Anchor {
                                    query_pos: seed_pos,
                                    ref_pos: pos,
                                    len: self.k,
                                }]
                            });
                    }
                }
            }
        };

        seed(&query, Strand::Forward);
        seed(&rc_query, Strand::Reverse);

        ref_to_anchors
    }

    fn calc_best_mapping<'b>(
        &mut self,
        seq_id: SequenceId,
        strand: Strand,
        query: &[u8],
        rev_query: &[u8],
        chains: impl Iterator<Item = &'b Chain>,
        score_threshold: i32,
    ) -> Option<Mapping> {
        let seq_range = self.index.seq_range(seq_id);
        let seq = &self.index.seq[seq_range.clone()];
        let rev_seq = sequence::reverse(seq);

        let align_config = self.aligner.config();
        let match_score = align_config.match_score as i32;
        let gap_open_penalty = align_config.gap_open_penalty as i32;
        let gap_extend_penalty = align_config.gap_extend_penalty as i32;

        let calc_bandwidth = |remaining_len, score, best_score| -> i32 {
            let max_gaps: i32 =
                (score + match_score * remaining_len as i32 - best_score - gap_open_penalty)
                    / gap_extend_penalty;
            (max_gaps + 1).max(1) + 1
        };

        let mut best_mapping = None;
        let mut best_score = score_threshold;

        'chain_loop: for chain in chains {
            let anchors = &chain.anchors;
            let mut score = 0;
            let mut remaining_len = query.len();

            // exact matches (anchors)
            let sum_anchors_len: usize = anchors.iter().map(|anchor| anchor.len).sum();
            score += match_score * sum_anchors_len as i32;
            remaining_len -= sum_anchors_len;

            // between anchors
            for (a, b) in anchors.iter().tuple_windows() {
                score += self.aligner.banded_global_align(
                    &query[a.query_pos + a.len..b.query_pos],
                    &self.index.seq[a.ref_pos + a.len..b.ref_pos],
                    calc_bandwidth(remaining_len, score, best_score),
                );

                remaining_len -= b.query_pos - a.query_pos - a.len;
                if score + match_score * remaining_len as i32 <= best_score {
                    continue 'chain_loop;
                }
            }

            let first_anchor = &anchors[0];
            let last_anchor = &anchors[anchors.len() - 1];

            let query_anchor_start = first_anchor.query_pos;
            let query_anchor_end = last_anchor.query_pos + last_anchor.len;

            let ref_anchor_start = first_anchor.ref_pos;
            let ref_anchor_end = last_anchor.ref_pos + last_anchor.len;

            // left extension
            if query_anchor_start > 0 {
                let bandwidth = calc_bandwidth(remaining_len, score, best_score);
                let target = if ref_anchor_start > seq_range.start {
                    &rev_seq[rev_seq.len() - (ref_anchor_start - seq_range.start) + 1..]
                } else {
                    &[]
                };

                score += self.aligner.banded_extension_align(
                    &rev_query[rev_query.len() - query_anchor_start + 1..],
                    &target[..target.len().min(query_anchor_start + bandwidth as usize)],
                    bandwidth,
                );

                remaining_len -= query_anchor_start;
                if score + match_score * remaining_len as i32 <= best_score {
                    continue;
                }
            }

            // right extension
            if query_anchor_end < query.len() {
                let bandwidth = calc_bandwidth(remaining_len, score, best_score);
                let target = if ref_anchor_end < seq_range.end {
                    &seq[ref_anchor_end - seq_range.start..]
                } else {
                    &[]
                };

                score += self.aligner.banded_extension_align(
                    &query[query_anchor_end..],
                    &target[..target
                        .len()
                        .min(query.len() - query_anchor_end + bandwidth as usize)],
                    bandwidth,
                );

                remaining_len -= query.len() - query_anchor_end;
            }

            assert_eq!(remaining_len, 0);

            if score > best_score {
                best_mapping = Some(Mapping {
                    seq_id,
                    pos: (ref_anchor_start - seq_range.start).saturating_sub(query_anchor_start),
                    strand,
                    score,
                });
                best_score = score;
            }
        }

        best_mapping
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
    use crate::{index::IndexBuilder, sequence};
    use itertools::Itertools;

    fn check_map(seqs: &[(&[u8], &[u8])], query: &[u8], expected: &[(&[u8], Strand)]) {
        let mut fasta: Vec<u8> = Vec::new();
        for (name, seq) in seqs {
            fasta.push(b'>');
            fasta.extend(*name);
            fasta.push(b'\n');
            fasta.extend(*seq);
            fasta.push(b'\n');
        }

        let cursor = std::io::Cursor::new(&fasta);
        let index = IndexBuilder::new(cursor).build().unwrap();
        let mut mapper = MapperBuilder::new(&index)
            .k(3)
            .min_score_fraction(0.25)
            .build();

        let query = sequence::encode(&query);

        let got: Vec<_> = mapper
            .map(&query)
            .iter()
            .map(|mapping| (mapping.seq_id, mapping.strand))
            .unique()
            .map(|(seq_id, strand)| (index.seq_name(seq_id), strand))
            .sorted()
            .collect();
        let expected: Vec<_> = expected.iter().sorted().copied().collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn map() {
        check_map(
            &[(b"foo", b"agctagt")],
            b"gct",
            &[(b"foo", Strand::Forward), (b"foo", Strand::Reverse)],
        );
        check_map(
            &[(b"foo", b"agctagt")],
            b"gctagt",
            &[(b"foo", Strand::Forward)],
        );
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"ggatcgatggagctctt",
            &[(b"foo", Strand::Forward)],
        );
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"tcgatataggagcgctta",
            &[(b"foo", Strand::Forward)],
        );
        check_map(
            &[(b"foo", b"agct"), (b"bar", b"tgta")],
            b"gct",
            &[(b"foo", Strand::Forward), (b"foo", Strand::Reverse)],
        );
        check_map(
            &[
                (b"foo", b"atcgggatatatggagagcttagag"),
                (b"bar", b"tatggagggagccttag"),
            ],
            b"atggagagcttag",
            &[(b"foo", Strand::Forward), (b"bar", Strand::Forward)],
        );

        check_map(
            &[(b"foo", b"agctagt")],
            b"agc",
            &[(b"foo", Strand::Reverse), (b"foo", Strand::Forward)],
        );
        check_map(
            &[(b"foo", b"agctagt")],
            b"actagc",
            &[(b"foo", Strand::Reverse)],
        );
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"aagagctccatcgatcc",
            &[(b"foo", Strand::Reverse)],
        );
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"taagcgctcctatatcga",
            &[(b"foo", Strand::Reverse)],
        );
        check_map(
            &[(b"foo", b"agct"), (b"bar", b"tgta")],
            b"gct",
            &[(b"foo", Strand::Forward), (b"foo", Strand::Reverse)],
        );
        check_map(
            &[
                (b"foo", b"atcgggatatatggagagcttagag"),
                (b"bar", b"tatggagggagccttag"),
            ],
            b"ctaagctctccat",
            &[(b"foo", Strand::Reverse), (b"bar", Strand::Reverse)],
        );
    }
}
