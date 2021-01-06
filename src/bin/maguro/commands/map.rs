mod parallel;

use super::Command;
use anyhow::{anyhow, Result};
use bio::io::fasta::{self, FastaRead};
use maguro::{
    index::Index,
    mapper::{align::AlignmentConfig, LibraryType, Mapper, MapperBuilder},
    sam, sequence, utils,
};
use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Write},
    path::PathBuf,
};
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
pub struct MapCommand {
    #[structopt(short, long)]
    index: PathBuf,

    #[structopt(short, long, conflicts_with_all(&["mate1", "mate2"]))]
    reads: Option<PathBuf>,
    #[structopt(short = "1", long)]
    mates1: Option<PathBuf>,
    #[structopt(short = "2", long)]
    mates2: Option<PathBuf>,

    #[structopt(short, long, default_value = "fr-unstranded")]
    library_type: LibraryType,
    #[structopt(long, default_value = "1000")]
    max_fragment_len: usize,

    #[structopt(short = "k", long, default_value = "31")]
    seed_min_len: usize,
    #[structopt(short, default_value = "10")]
    multiplicity: usize,

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
    #[structopt(short, long, default_value = "1")]
    threads: usize,
}

impl Command for MapCommand {
    fn run(self) -> Result<()> {
        let index: Index = {
            let reader = BufReader::new(File::open(&self.index)?);
            bincode::deserialize_from(reader)?
        };

        let mut out = io::stdout();
        {
            let mut out = BufWriter::new(out.lock());
            sam::write_header(&mut out, &index)?;
            out.flush()?;
        }

        if self.threads > 1 {
            parallel::main(self, &mut out, &index)?;
            return Ok(());
        }

        let mut mapper = build_mapper(&index, &self);
        let mut out = BufWriter::new(out.lock());

        if let Some(ref reads) = self.reads {
            // single-end

            let mut reader = fasta::Reader::from_file(&reads)?;
            let mut record = fasta::Record::new();
            reader.read(&mut record)?;

            while !record.is_empty() {
                record.check().map_err(|e| anyhow!(e.to_owned()))?;

                let qname = utils::extract_name_bytes(record.id(), &self.header_sep);
                map_single(&mut out, &index, &mut mapper, qname, record.seq())?;

                reader.read(&mut record)?;
            }
        } else {
            // paired-end

            let mut reader1 = fasta::Reader::from_file(self.mates1.as_ref().unwrap())?;
            let mut reader2 = fasta::Reader::from_file(self.mates2.as_ref().unwrap())?;

            let mut record1 = fasta::Record::new();
            let mut record2 = fasta::Record::new();
            reader1.read(&mut record1)?;
            reader2.read(&mut record2)?;

            while !record1.is_empty() {
                record1.check().map_err(|e| anyhow!(e.to_owned()))?;
                record2.check().map_err(|e| anyhow!(e.to_owned()))?;

                let qname = utils::extract_name_bytes(record1.id(), &self.header_sep);
                map_pair(
                    &mut out,
                    &index,
                    &mut mapper,
                    qname,
                    record1.seq(),
                    record2.seq(),
                )?;

                reader1.read(&mut record1)?;
                reader2.read(&mut record2)?;
            }

            assert!(record2.is_empty());
        }

        out.flush()?;
        Ok(())
    }
}

fn build_mapper<'a>(index: &'a Index, config: &MapCommand) -> Mapper<'a> {
    MapperBuilder::new(&index)
        .library_type(config.library_type)
        .seed_min_len(config.seed_min_len)
        .seed_max_hits(config.multiplicity)
        .consensus_fraction(config.consensus_fraction)
        .coverage_score_ratio(config.coverage_score_ratio)
        .max_splice_gap(config.max_splice_gap)
        .min_score_fraction(config.min_score_fraction)
        .alignment_config(config.alignment_config.clone())
        .build()
}

fn map_single<'a, W: Write>(
    mut out: W,
    index: &Index,
    mapper: &mut Mapper<'a>,
    qname: &[u8],
    seq: &[u8],
) -> Result<()> {
    let encoded_seq = sequence::encode(&seq);

    let mut mappings = mapper.map_single(&encoded_seq);
    if mappings.is_empty() {
        sam::write_unmapped_single(&mut out, qname, seq)?;
    } else {
        mappings.sort_by(|a, b| b.score.cmp(&a.score));

        let mut secondary = false;
        let mut rc_query_cache = None;

        for mapping in mappings {
            sam::write_mapping_single(
                &mut out,
                &index,
                qname,
                seq,
                &mapping,
                secondary,
                &mut rc_query_cache,
            )?;
            secondary = true;
        }
    }

    Ok(())
}

fn map_pair<'a, W: Write>(
    mut out: W,
    index: &Index,
    mapper: &mut Mapper<'a>,
    qname: &[u8],
    seq1: &[u8],
    seq2: &[u8],
) -> Result<()> {
    let encoded_seq1 = sequence::encode(&seq1);
    let encoded_seq2 = sequence::encode(&seq2);

    let mut mappings = mapper.map_pair(&encoded_seq1, &encoded_seq2);
    if mappings.is_empty() {
        sam::write_unmapped_pair(&mut out, qname, seq1, seq2)?;
    } else {
        mappings.sort_by(|a, b| b.score.cmp(&a.score));

        let mut secondary = false;
        let mut rc_query_cache1 = None;
        let mut rc_query_cache2 = None;

        for mapping in mappings {
            sam::write_mapping_pair(
                &mut out,
                &index,
                qname,
                seq1,
                seq2,
                &mapping,
                secondary,
                &mut rc_query_cache1,
                &mut rc_query_cache2,
            )?;
            secondary = true;
        }
    }

    Ok(())
}
