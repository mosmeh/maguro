use super::{Anchor, Mapper, Strand};
use crate::{index::SequenceId, sequence};
use std::collections::HashMap;

pub struct PairMapping {
    pub seq_id: SequenceId,
    pub pos1: usize,
    pub pos2: usize,
    pub strand: Strand,
    pub score: i32,
}

impl Mapper<'_> {
    pub fn map_pair(&mut self, query1: &[u8], query2: &[u8]) -> Vec<PairMapping> {
        if query1.len() < self.seed_min_len || query2.len() < self.seed_min_len {
            // TODO
            return Vec::new();
        }

        // seeding
        let rc_query1 = sequence::reverse_complement(&query1);
        let ref_to_anchors1 = self.search_anchors(query1, &rc_query1);
        if ref_to_anchors1.is_empty() {
            return Vec::new();
        }

        let rc_query2 = sequence::reverse_complement(&query2);
        let ref_to_anchors2 = self.search_anchors(query2, &rc_query2);
        if ref_to_anchors2.is_empty() {
            return Vec::new();
        }

        // chaining
        let ref_to_chains1 = self.chain_anchors(ref_to_anchors1);
        if ref_to_chains1.is_empty() {
            return Vec::new();
        }

        let ref_to_chains2 = self.chain_anchors(ref_to_anchors2);
        if ref_to_chains2.is_empty() {
            return Vec::new();
        }

        // match pairs
        let mut best_pair_score = f64::MIN;
        let mut pair_score_threshold = f64::MIN;
        let mut ref_to_pairs = HashMap::new();

        // "IU" in https://salmon.readthedocs.io/en/latest/library_type.html
        for ((seq_id, strand), chains1) in &ref_to_chains1 {
            if let Some(chains2) = &ref_to_chains2.get(&(*seq_id, strand.opposite())) {
                let mut pairs = Vec::new();

                for chain1 in chains1 {
                    let first_anchor1 = &chain1.anchors[0];
                    let last_anchor1 = &chain1.anchors[chain1.anchors.len() - 1];

                    for chain2 in *chains2 {
                        let first_anchor2 = &chain2.anchors[0];
                        let last_anchor2 = &chain2.anchors[chain2.anchors.len() - 1];

                        let fragment_len = if first_anchor1.ref_pos < first_anchor2.ref_pos {
                            last_anchor2.ref_pos + last_anchor2.len - first_anchor1.ref_pos
                        } else {
                            last_anchor1.ref_pos + last_anchor1.len - first_anchor2.ref_pos
                        };
                        if fragment_len > self.max_fragment_len {
                            continue;
                        }

                        let score = chain1.score + chain2.score;
                        if score > best_pair_score {
                            best_pair_score = score;
                            pair_score_threshold = self.coverage_score_ratio * best_pair_score;
                        }

                        if score >= pair_score_threshold {
                            pairs.push((chain1, chain2));
                        }
                    }
                }

                ref_to_pairs.insert((seq_id, strand), pairs);
            }
        }

        // base-to-base alignment
        let rev_query1 = sequence::reverse(query1);
        let compl_query1 = sequence::complement(query1);

        let rev_query2 = sequence::reverse(query2);
        let compl_query2 = sequence::complement(query2);

        let match_score = self.aligner.config().match_score as i32;
        let max_align_score1 = query1.len() as i32 * match_score;
        let align_score_threshold1 = (max_align_score1 as f64 * self.min_score_fraction) as i32;
        let max_align_score2 = query2.len() as i32 * match_score;
        let align_score_threshold2 = (max_align_score2 as f64 * self.min_score_fraction) as i32;

        let mut anchors_to_scores1: HashMap<(SequenceId, Strand, &[Anchor]), i32> = HashMap::new();
        let mut anchors_to_scores2: HashMap<(SequenceId, Strand, &[Anchor]), i32> = HashMap::new();

        for ((seq_id, strand), pairs) in ref_to_pairs.iter_mut() {
            pairs.retain(|pair| pair.0.score + pair.1.score >= pair_score_threshold);

            let seq_range = self.index.seq_range(**seq_id);
            let seq = &self.index.seq[seq_range.clone()];
            let rev_seq = sequence::reverse(seq);

            let (align_query1, rev_align_query1): (&[u8], &[u8]) = match strand {
                Strand::Forward => (query1, &rev_query1),
                Strand::Reverse => (&rc_query1, &compl_query1),
            };
            let (align_query2, rev_align_query2): (&[u8], &[u8]) = match strand.opposite() {
                Strand::Forward => (query2, &rev_query2),
                Strand::Reverse => (&rc_query2, &compl_query2),
            };

            for pair in pairs {
                anchors_to_scores1
                    .entry((**seq_id, **strand, &pair.0.anchors))
                    .or_insert_with(|| {
                        self.calc_align_score(
                            align_query1,
                            rev_align_query1,
                            seq,
                            &rev_seq,
                            &seq_range,
                            pair.0,
                            align_score_threshold1,
                        )
                        .unwrap_or(i32::MIN)
                    });
                anchors_to_scores2
                    .entry((**seq_id, **strand, &pair.1.anchors))
                    .or_insert_with(|| {
                        self.calc_align_score(
                            align_query2,
                            rev_align_query2,
                            seq,
                            &rev_seq,
                            &seq_range,
                            pair.1,
                            align_score_threshold2,
                        )
                        .unwrap_or(i32::MIN)
                    });
            }
        }

        let mut mappings = Vec::new();

        for ((seq_id, strand), pairs) in ref_to_pairs {
            let seq_range = self.index.seq_range(*seq_id);

            let mut best_mapping = None;
            let mut best_score = i32::MIN;

            for pair in pairs {
                let score1 = anchors_to_scores1[&(*seq_id, *strand, pair.0.anchors.as_slice())];
                let score2 = anchors_to_scores2[&(*seq_id, *strand, pair.1.anchors.as_slice())];

                if score1 > i32::MIN && score2 > i32::MIN && score1 + score2 > best_score {
                    let score = score1 + score2;
                    best_score = score;

                    let first_anchor1 = &pair.0.anchors[0];
                    let first_anchor2 = &pair.1.anchors[0];

                    best_mapping = Some(PairMapping {
                        seq_id: *seq_id,
                        pos1: (first_anchor1.ref_pos - seq_range.start)
                            .saturating_sub(first_anchor1.query_pos),
                        pos2: (first_anchor2.ref_pos - seq_range.start)
                            .saturating_sub(first_anchor2.query_pos),
                        strand: *strand,
                        score,
                    });
                }
            }

            if let Some(mapping) = best_mapping {
                mappings.push(mapping);
            }
        }

        mappings
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        index::IndexBuilder,
        mapper::{MapperBuilder, Strand},
        sequence,
    };
    use itertools::Itertools;

    fn check_map(
        seqs: &[(&[u8], &[u8])],
        query1: &[u8],
        query2: &[u8],
        expected: &[(&[u8], Strand)],
    ) {
        let mut fasta: Vec<u8> = Vec::new();
        for (name, seq) in seqs {
            fasta.push(b'>');
            fasta.extend(*name);
            fasta.push(b'\n');
            fasta.extend(*seq);
            fasta.push(b'\n');
        }

        let cursor = std::io::Cursor::new(&fasta);
        let index = IndexBuilder::new(cursor).bucket_width(2).build().unwrap();
        let mut mapper = MapperBuilder::new(&index)
            .max_fragment_len(15)
            .seed_min_len(3)
            .min_score_fraction(0.25)
            .build();

        let query1 = sequence::encode(&query1);
        let query2 = sequence::encode(&query2);

        let got: Vec<_> = mapper
            .map_pair(&query1, &query2)
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
            b"act",
            &[(b"foo", Strand::Forward)],
        );
        check_map(&[(b"foo", b"agctagt")], b"gctagt", b"gctagt", &[]);
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"ggatcgatggagctctt",
            b"ctctaagc",
            &[],
        );
        check_map(
            &[(b"foo", b"atcgggatatatggagagcttagag")],
            b"ggatcgatggagctctt",
            b"catata",
            &[(b"foo", Strand::Forward)],
        );
        check_map(
            &[
                (b"foo", b"atcgggatatatggagagcttagag"),
                (b"bar", b"tatggagggagccttag"),
            ],
            b"ctaagctctccat",
            b"tatgg",
            &[(b"foo", Strand::Reverse)],
        );
    }
}
