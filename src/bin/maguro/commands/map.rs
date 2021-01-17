mod parallel;
mod serial;

use super::Command;
use anyhow::Result;
use maguro::{
    index::Index,
    mapper::{align::AlignmentConfig, LibraryType, Mapper, MapperBuilder},
    sequence,
};
use std::{
    fs::File,
    io::{BufReader, Write},
    path::PathBuf,
    time::Instant,
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

    #[structopt(long, default_value = env!("CARGO_PKG_NAME"))]
    sam_pg: String,

    #[structopt(short, long, default_value = "1")]
    threads: usize,
    #[structopt(short, long, default_value = "1")]
    chunk: usize,
}

impl Command for MapCommand {
    fn run(self) -> Result<()> {
        eprintln!("{:#?}", self);

        eprintln!("Loading index");
        let index: Index = {
            let reader = BufReader::new(File::open(&self.index)?);
            bincode::deserialize_from(reader)?
        };

        let mapper = MapperBuilder::new(&index)
            .library_type(self.library_type)
            .seed_min_len(self.seed_min_len)
            .seed_max_hits(self.multiplicity)
            .consensus_fraction(self.consensus_fraction)
            .coverage_score_ratio(self.coverage_score_ratio)
            .max_splice_gap(self.max_splice_gap)
            .min_score_fraction(self.min_score_fraction)
            .alignment_config(self.alignment_config.clone())
            .build();

        let start_time = Instant::now();

        if self.threads > 1 {
            parallel::main(self, &index, &mapper)?;
        } else {
            serial::main(self, &index, &mapper)?;
        }

        eprintln!("Elapsed {}ms", start_time.elapsed().as_millis());
        eprintln!("Finished");
        eprintln!(
            "#do_one_step = {}",
            maguro::index::suffix_array::STEP_COUNTER.load(std::sync::atomic::Ordering::Relaxed)
        );
        let found =
            maguro::index::suffix_array::FOUND_COUNTER.load(std::sync::atomic::Ordering::Relaxed);
        let search =
            maguro::index::suffix_array::SEARCH_COUNTER.load(std::sync::atomic::Ordering::Relaxed);
        eprintln!(
            "Found {} / {} ({:.2}%)",
            found,
            search,
            found as f64 * 100.0 / search as f64
        );

        Ok(())
    }
}

fn map_single<'a, W: Write>(
    _out: W,
    _index: &Index,
    mapper: &Mapper<'a>,
    _qname: &[u8],
    seq: &[u8],
) -> Result<bool> {
    let encoded_seq = sequence::encode(&seq);

    let mappings = mapper.map_single(&encoded_seq);
    Ok(!mappings.is_empty())
}

fn map_pair<'a, W: Write>(
    _out: W,
    _index: &Index,
    mapper: &Mapper<'a>,
    _qname: &[u8],
    seq1: &[u8],
    seq2: &[u8],
) -> Result<bool> {
    let encoded_seq1 = sequence::encode(&seq1);
    let encoded_seq2 = sequence::encode(&seq2);

    let mappings = mapper.map_pair(&encoded_seq1, &encoded_seq2);
    Ok(!mappings.is_empty())
}
