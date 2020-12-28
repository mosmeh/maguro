use super::Command;
use bio::io::fasta::{self, FastaRead};
use maguro::{
    index::Index,
    mapper::{align::AlignmentConfig, MapperBuilder},
    sam::SamWriter,
    utils,
};
use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::PathBuf,
};
use structopt::StructOpt;

#[derive(StructOpt)]
pub struct MapCommand {
    #[structopt(short, long)]
    index: PathBuf,
    #[structopt(short, long)]
    read: PathBuf,
    #[structopt(short, default_value = "31")]
    k: usize,
    #[structopt(long, default_value = "0.65")]
    consensus_fraction: f64,
    #[structopt(long, default_value = "0.6")]
    coverage_score_ratio: f64,
    #[structopt(long, default_value = "100")]
    max_splice_gap: usize,
    #[structopt(long, default_value = "0.65")]
    min_score_fraction: f64,
    #[structopt(flatten)]
    alignment_config: AlignmentConfig,
    #[structopt(long)]
    header_sep: Option<String>,
}

impl Command for MapCommand {
    fn run(self) -> anyhow::Result<()> {
        let index: Index = {
            let reader = BufReader::new(File::open(self.index)?);
            bincode::deserialize_from(reader)?
        };
        let mapper = MapperBuilder::new(&index)
            .k(self.k)
            .consensus_fraction(self.consensus_fraction)
            .coverage_score_ratio(self.coverage_score_ratio)
            .max_splice_gap(self.max_splice_gap)
            .min_score_fraction(self.min_score_fraction)
            .alignment_config(self.alignment_config)
            .build();

        let out = std::io::stdout();
        let out = BufWriter::new(out.lock());
        let mut sam_writer = SamWriter::new(out);
        sam_writer.write_header(&index)?;

        let mut reader = fasta::Reader::from_file(self.read).expect("Failed to open FASTA file");
        let mut record = fasta::Record::new();
        reader.read(&mut record).expect("Failed to parse record");

        while !record.is_empty() {
            record.check().expect("Invalid record");

            let qname = utils::extract_byte_name(record.id(), &self.header_sep);
            let query = utils::encode_seq(&record.seq());

            let mappings = mapper.map(&query);
            if mappings.is_empty() {
                sam_writer.write_unmapped(qname, record.seq())?;
            } else {
                let mut secondary = false;
                for mapping in mappings {
                    sam_writer.write_mapping(&index, qname, &record.seq(), &mapping, secondary)?;
                    secondary = true;
                }
            }

            reader.read(&mut record).expect("Failed to parse record");
        }

        Ok(())
    }
}