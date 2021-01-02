mod parallel;

use super::Command;
use anyhow::{anyhow, Result};
use bio::io::fasta::{self, FastaRead};
use maguro::{
    index::Index,
    mapper::{align::AlignmentConfig, Mapper, MapperBuilder},
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
    #[structopt(short, long)]
    read: PathBuf,
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

        let mut reader = fasta::Reader::from_file(&self.read)?;

        if self.threads > 1 {
            parallel::parallel_map(self, &mut out, &index, &mut reader)?;
            return Ok(());
        }

        let mut mapper = build_mapper(&index, &self);
        let mut out = BufWriter::new(out.lock());

        let mut record = fasta::Record::new();
        reader.read(&mut record)?;

        while !record.is_empty() {
            record.check().map_err(|e| anyhow!(e.to_owned()))?;

            let qname = utils::extract_name_bytes(record.id(), &self.header_sep);
            map(&mut out, &index, &mut mapper, qname, record.seq())?;

            reader.read(&mut record)?;
        }

        out.flush()?;
        Ok(())
    }
}

fn build_mapper<'a>(index: &'a Index, config: &MapCommand) -> Mapper<'a> {
    MapperBuilder::new(&index)
        .seed_min_len(config.seed_min_len)
        .seed_max_hits(config.multiplicity)
        .consensus_fraction(config.consensus_fraction)
        .coverage_score_ratio(config.coverage_score_ratio)
        .max_splice_gap(config.max_splice_gap)
        .min_score_fraction(config.min_score_fraction)
        .alignment_config(config.alignment_config.clone())
        .build()
}

fn map<'a, W: Write>(
    mut out: W,
    index: &Index,
    mapper: &mut Mapper<'a>,
    qname: &[u8],
    seq: &[u8],
) -> Result<()> {
    let encoded_seq = sequence::encode(&seq);

    let mut mappings = mapper.map(&encoded_seq);
    if mappings.is_empty() {
        sam::write_unmapped(&mut out, qname, seq)?;
    } else {
        mappings.sort_by(|a, b| b.score.cmp(&a.score));

        let mut secondary = false;
        let mut rc_query_cache = None;

        for mapping in mappings {
            sam::write_mapping(
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
