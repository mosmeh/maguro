pub mod align;
mod chain;
pub mod pair;
pub mod single;

use crate::index::{Index, SequenceId};
use align::{Aligner, AlignmentConfig};
use chain::Chain;
use itertools::Itertools;
use rustc_hash::FxHashMap;
use std::ops::Range;

#[derive(Debug, Clone, Copy)]
pub enum LibraryType {
    Unstranded,
    FirstStrand,
    SecondStrand,
}

impl std::str::FromStr for LibraryType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        match &s.to_lowercase()[..] {
            "fr-unstranded" => Ok(Self::Unstranded),
            "fr-firststrand" => Ok(Self::FirstStrand),
            "fr-secondstrand" => Ok(Self::SecondStrand),
            _ => Err(format!(
                "Unknown library type {}. \
            Valid values are: fr-unstranded, fr-firststrand, fr-secondstrand",
                s
            )),
        }
    }
}

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

    pub fn opposite(self) -> Strand {
        match self {
            Self::Forward => Self::Reverse,
            Self::Reverse => Self::Forward,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Anchor {
    query_pos: usize,
    ref_pos: usize,
    len: usize,
}

pub struct Mapper<'a> {
    index: &'a Index,
    library_type: LibraryType,
    seed_min_len: usize,
    seed_max_hits: usize,
    max_fragment_len: usize,
    consensus_fraction: f64,
    coverage_score_ratio: f64,
    max_splice_gap: usize,
    min_score_fraction: f64,
    aligner: Aligner,
}

impl Mapper<'_> {
    fn search_anchors(
        &self,
        query: &[u8],
        rc_query: &[u8],
        is_read1: bool,
    ) -> FxHashMap<(SequenceId, Strand), Vec<Anchor>> {
        let mut ref_to_anchors: FxHashMap<(SequenceId, Strand), Vec<Anchor>> = FxHashMap::default();

        let mut seed = |query: &[u8], strand| {
            for seed_pos in 0..=(query.len() - self.seed_min_len) {
                let result = self.index.sa.extension_search(
                    &self.index.seq,
                    &query[seed_pos..],
                    self.seed_min_len,
                    self.seed_max_hits,
                );
                if let Some((range, len)) = result {
                    for i in range {
                        let pos = self.index.sa.array[i] as usize;
                        let id = self.index.seq_id_from_pos(pos);

                        ref_to_anchors
                            .entry((id, strand))
                            .and_modify(|anchors| {
                                anchors.push(Anchor {
                                    query_pos: seed_pos,
                                    ref_pos: pos,
                                    len,
                                })
                            })
                            .or_insert_with(|| {
                                vec![Anchor {
                                    query_pos: seed_pos,
                                    ref_pos: pos,
                                    len,
                                }]
                            });
                    }
                }
            }
        };

        match (self.library_type, is_read1) {
            (LibraryType::Unstranded, _) => {
                seed(&query, Strand::Forward);
                seed(&rc_query, Strand::Reverse);
            }
            (LibraryType::FirstStrand, false) | (LibraryType::SecondStrand, true) => {
                seed(&query, Strand::Forward);
            }
            (LibraryType::SecondStrand, false) | (LibraryType::FirstStrand, true) => {
                seed(&rc_query, Strand::Reverse);
            }
        }

        ref_to_anchors
    }

    fn chain_anchors(
        &self,
        ref_to_anchors: FxHashMap<(SequenceId, Strand), Vec<Anchor>>,
    ) -> FxHashMap<(SequenceId, Strand), Vec<Chain>> {
        let mut ref_to_chains = FxHashMap::default();
        let mut best_score = f64::MIN;
        let mut score_threshold = f64::MIN;

        for ((seq_id, strand), anchors) in ref_to_anchors {
            let mut chains = chain::chain_anchors(anchors, self.max_splice_gap, score_threshold);
            if chains.is_empty() {
                continue;
            }

            let ref_best_score = chains
                .iter()
                .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
                .unwrap()
                .score;
            if ref_best_score > best_score {
                best_score = ref_best_score;
                score_threshold = best_score * self.consensus_fraction;
                chains.retain(|chain| chain.score >= score_threshold);
            }

            if !chains.is_empty() {
                ref_to_chains.insert((seq_id, strand), chains);
            }
        }

        for (_, chains) in ref_to_chains.iter_mut() {
            chains.retain(|chain| chain.score >= score_threshold);
        }

        ref_to_chains
    }

    #[allow(clippy::too_many_arguments)]
    fn calc_align_score(
        &self,
        query: &[u8],
        rev_query: &[u8],
        seq: &[u8],
        rev_seq: &[u8],
        seq_range: &Range<usize>,
        chain: &Chain,
        score_threshold: i32,
    ) -> Option<i32> {
        let align_config = self.aligner.config();
        let match_score = align_config.match_score as i32;

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
                self.calc_bandwidth(remaining_len, score, score_threshold),
            );

            remaining_len -= b.query_pos - a.query_pos - a.len;
            if (score + match_score * remaining_len as i32) < score_threshold {
                return None;
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
            let bandwidth = self.calc_bandwidth(remaining_len, score, score_threshold);
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
            if (score + match_score * remaining_len as i32) < score_threshold {
                return None;
            }
        }

        // right extension
        if query_anchor_end < query.len() {
            let bandwidth = self.calc_bandwidth(remaining_len, score, score_threshold);
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

        if score >= score_threshold {
            Some(score)
        } else {
            None
        }
    }

    fn calc_bandwidth(&self, remaining_len: usize, score: i32, score_threshold: i32) -> i32 {
        let align_config = self.aligner.config();
        let match_score = align_config.match_score as i32;
        let gap_open_penalty = align_config.gap_open_penalty as i32;
        let gap_extend_penalty = align_config.gap_extend_penalty as i32;

        let max_gaps: i32 =
            (score + match_score * remaining_len as i32 - score_threshold - gap_open_penalty)
                / gap_extend_penalty;
        (max_gaps + 1).max(1) + 1
    }
}

pub struct MapperBuilder<'a> {
    index: &'a Index,
    library_type: LibraryType,
    seed_min_len: usize,
    seed_max_hits: usize,
    max_fragment_len: usize,
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
            library_type: LibraryType::Unstranded,
            seed_min_len: 31,
            seed_max_hits: 10,
            max_fragment_len: 1000,
            consensus_fraction: 0.65,
            coverage_score_ratio: 0.6,
            max_splice_gap: 100,
            min_score_fraction: 0.65,
            alignment_config: Default::default(),
        }
    }

    pub fn library_type(&mut self, library_type: LibraryType) -> &mut Self {
        self.library_type = library_type;
        self
    }

    pub fn seed_min_len(&mut self, seed_min_len: usize) -> &mut Self {
        self.seed_min_len = seed_min_len;
        self
    }

    pub fn seed_max_hits(&mut self, seed_max_hits: usize) -> &mut Self {
        self.seed_max_hits = seed_max_hits;
        self
    }

    pub fn max_fragment_len(&mut self, max_fragment_len: usize) -> &mut Self {
        self.max_fragment_len = max_fragment_len;
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
            library_type: self.library_type,
            seed_min_len: self.seed_min_len,
            seed_max_hits: self.seed_max_hits,
            max_fragment_len: self.max_fragment_len,
            consensus_fraction: self.consensus_fraction,
            coverage_score_ratio: self.coverage_score_ratio,
            max_splice_gap: self.max_splice_gap,
            min_score_fraction: self.min_score_fraction,
            aligner: Aligner::new(self.alignment_config.clone()),
        }
    }
}
