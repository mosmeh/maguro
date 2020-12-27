mod chain;

use crate::index::{Index, SequenceId};
use chain::Chain;
use std::collections::HashMap;

#[derive(Clone)]
pub struct Anchor {
    query_pos: usize,
    ref_pos: usize,
    len: usize,
}

pub struct Mapper<'a> {
    index: &'a Index,
    k: usize,
    consensus_fraction: f64,
    coverage_score_ratio: f64,
    max_splice_gap: usize,
}

impl Mapper<'_> {
    pub fn map(&self, query: &[u8]) -> Vec<(SequenceId, usize)> {
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
            let mut chains = chain::chain_anchors(anchors, self.max_splice_gap);
            let best_score = chains
                .iter()
                .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
                .unwrap()
                .score;
            global_best_score = global_best_score.max(best_score);

            let min_required_score = best_score * self.coverage_score_ratio;
            chains.retain(|chain| chain.score >= min_required_score);

            ref_chains.insert(seq_id, chains);
        }

        let min_required_score = global_best_score * self.consensus_fraction;
        for (_, chains) in ref_chains.iter_mut() {
            chains.retain(|chain| chain.score >= min_required_score);
        }

        // TODO: base-to-base alignments

        ref_chains
            .into_iter()
            .map(|(seq_id, chains)| {
                let best_chain = chains
                    .into_iter()
                    .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
                    .unwrap();
                let pos = best_chain
                    .anchors
                    .iter()
                    .min_by(|a, b| a.ref_pos.cmp(&b.ref_pos))
                    .unwrap()
                    .ref_pos;
                (seq_id, pos - self.index.seq_range(seq_id).start)
            })
            .collect()
    }
}

pub struct MapperBuilder<'a> {
    index: &'a Index,
    k: usize,
    consensus_fraction: f64,
    coverage_score_ratio: f64,
    max_splice_gap: usize,
}

impl<'a> MapperBuilder<'a> {
    pub fn new(index: &'a Index) -> Self {
        Self {
            index,
            k: 15,
            consensus_fraction: 0.65,
            coverage_score_ratio: 0.6,
            max_splice_gap: 100,
        }
    }

    pub fn k(mut self, k: usize) -> Self {
        self.k = k;
        self
    }

    pub fn consensus_fraction(mut self, consensus_fraction: f64) -> Self {
        self.consensus_fraction = consensus_fraction;
        self
    }

    pub fn coverage_score_ratio(mut self, coverage_score_ratio: f64) -> Self {
        self.coverage_score_ratio = coverage_score_ratio;
        self
    }

    pub fn max_splice_gap(mut self, max_splice_gap: usize) -> Self {
        self.max_splice_gap = max_splice_gap;
        self
    }

    pub fn build(self) -> Mapper<'a> {
        Mapper {
            index: self.index,
            k: self.k,
            consensus_fraction: self.consensus_fraction,
            coverage_score_ratio: self.coverage_score_ratio,
            max_splice_gap: self.max_splice_gap,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{index::IndexBuilder, utils};
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
        let mapper = MapperBuilder::new(&index).k(3).build();

        let mut query = query.to_owned();
        utils::encode_seq_in_place(&mut query);

        let got: Vec<_> = mapper
            .map(&query)
            .iter()
            .map(|(id, _)| index.seq_name(*id))
            .sorted()
            .collect();
        let expected: Vec<_> = expected.iter().sorted().copied().collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn map() {
        check_map(&[(b"foo", b"agctagt")], b"gct", &[b"foo"]);
        check_map(&[(b"foo", b"agctagt")], b"agctagta", &[b"foo"]);
        check_map(&[(b"foo", b"agctagt")], b"gca", &[]);
        check_map(&[(b"foo", b"atcgggatatatggagagcttagag")], b"gag", &[b"foo"]);
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"ataagagccct",
            &[b"foo"],
        );
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"ggatattggagagc",
            &[b"foo"],
        );
        check_map(&[(b"foo", b"agct"), (b"bar", b"tgta")], b"gct", &[b"foo"]);
        check_map(
            &[
                (b"foo", b"atcgggatatatggagagcttagag"),
                (b"bar", b"tttagagggttcttcgggatt"),
            ],
            b"gag",
            &[b"foo", b"bar"],
        );
    }
}
