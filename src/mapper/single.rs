use super::{Mapper, Strand};
use crate::{index::SequenceId, sequence};

pub struct SingleMapping {
    pub seq_id: SequenceId,
    pub pos: usize,
    pub strand: Strand,
    pub score: i32,
}

impl Mapper<'_> {
    pub fn map_single(&self, query: &[u8]) -> Vec<SingleMapping> {
        if query.len() < self.seed_min_len {
            // TODO
            return Vec::new();
        }

        let rc_query = sequence::reverse_complement(&query);

        // seeding
        let ref_to_anchors = self.search_anchors(query, &rc_query, true);
        if ref_to_anchors.is_empty() {
            return Vec::new();
        }

        // chaining
        let ref_to_chains = self.chain_anchors(ref_to_anchors);
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

        for ((seq_id, strand), chains) in ref_to_chains {
            let seq_range = self.index.seq_range(seq_id);
            let seq = &self.index.seq[seq_range.clone()];
            let rev_seq = sequence::reverse(seq);

            let (align_query, rev_align_query): (&[u8], &[u8]) = match strand {
                Strand::Forward => (query, &rev_query),
                Strand::Reverse => (&rc_query, &compl_query),
            };

            let mut best_mapping = None;
            let mut best_score = align_score_threshold;

            for chain in chains {
                let score = self.calc_align_score(
                    align_query,
                    &rev_align_query,
                    seq,
                    &rev_seq,
                    &seq_range,
                    &chain,
                    best_score,
                );
                if let Some(score) = score {
                    let first_anchor = &chain.anchors[0];
                    best_mapping = Some(SingleMapping {
                        seq_id,
                        pos: (first_anchor.ref_pos - seq_range.start)
                            .saturating_sub(first_anchor.query_pos),
                        strand,
                        score,
                    });
                    best_score = score;
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
        let index = IndexBuilder::new(cursor).bucket_width(2).build().unwrap();
        let mapper = MapperBuilder::new(&index)
            .seed_min_len(3)
            .min_score_fraction(0.25)
            .build();

        let query = sequence::encode(&query);

        let got: Vec<_> = mapper
            .map_single(&query)
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
