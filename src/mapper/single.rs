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

        let mut mappings = Vec::new();
        for ((seq_id, strand), anchors) in ref_to_anchors {
            mappings.push(SingleMapping {
                seq_id,
                pos: anchors[0].ref_pos,
                strand,
                score: 0,
            });
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
