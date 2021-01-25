use super::Command;
use maguro::index::{Index, SequenceId};
use prettytable::{format::consts::FORMAT_CLEAN, *};
use size::Size;
use std::{fs::File, io::BufReader, path::PathBuf};
use structopt::StructOpt;

#[derive(StructOpt)]
pub struct StatsCommand {
    #[structopt(short, long)]
    index: PathBuf,
}

impl Command for StatsCommand {
    #[allow(unused_assignments)]
    fn run(self) -> anyhow::Result<()> {
        let index: Index = {
            let reader = BufReader::new(File::open(&self.index)?);
            bincode::deserialize_from(reader)?
        };

        /*let num_seqs = index.num_seqs();
        let mut min_len = usize::MAX;
        let mut max_len = usize::MIN;
        let mut sum_len = 0;

        for i in 0..num_seqs {
            let len = index.seq(SequenceId(i)).len();
            min_len = min_len.min(len);
            max_len = max_len.max(len);
            sum_len += len;
        }

        println!("Sequence stats");
        let mut table = Table::new();
        table.set_format(*FORMAT_CLEAN);

        macro_rules! item {
            ($name:expr, $num:expr) => {
                table.add_row(row![$name, r -> $num]);
            }
        }
        item!("Num sequences", num_seqs);
        item!("Sum length", sum_len);
        item!("Min length", min_len);
        item!("Avg length", sum_len / num_seqs);
        item!("Max length", max_len);

        table.printstd();

        println!("\nIndex size");
        let mut table = Table::new();
        table.set_format(*FORMAT_CLEAN);

        let mut total_bytes = 0;
        macro_rules! from_bytes {
            ($name:expr, $bytes:expr) => {
                let formatted = format!("{}", Size::Bytes($bytes));
                table.add_row(row![$name, r -> formatted]);
                total_bytes += $bytes;
            }
        }
        macro_rules! from_slice {
            ($name:expr, $slice:expr) => {
                let bytes = $slice.len() * std::mem::size_of_val(&$slice[0]);
                from_bytes!($name, bytes);
            };
        }

        from_slice!("Sequences", &index.seq);
        from_slice!("Sequence boundaries", &index.ends);
        from_slice!("Sequence names", &index.name_arena);
        from_slice!("Sequence name boundaries", &index.name_ends);
        from_bytes!("Rank dictionary", index.rank_dict.size_bytes());
        table.add_empty_row();

        from_slice!("Suffix array", &index.sa.array);
        from_slice!("Child table", &index.sa.child);
        from_slice!("Buckets", &index.sa.offsets);
        table.add_empty_row();

        from_bytes!("Total", total_bytes);

        table.printstd();*/

        let mut non_zero = 0;
        for i in 0..(index.sa.offsets.len() - 1) {
            let len = index.sa.offsets[i + 1] - index.sa.offsets[i];
            if len > 0 {
                non_zero += 1;
            }
        }

        let mut kmers = std::collections::HashSet::new();
        for i in 0..=(index.seq.len() - index.sa.bucket_width) {
            let kmer = &index.seq[i..][..index.sa.bucket_width];
            if !kmer
                .iter()
                .any(|x| *x == 0 || *x == maguro::sequence::DUMMY_CODE)
            {
                kmers.insert(kmer);
            }
        }

        let mut counts = vec![0; index.sa.offsets.len() - 1];
        for kmer in &kmers {
            let mut idx = 0;
            for (j, x) in kmer.iter().enumerate() {
                idx |= (maguro::sequence::code_to_two_bit(*x) as usize) << (2 * j);
            }
            counts[idx] += 1;
        }
        let mut collision = 0;
        for count in counts {
            if count > 1 {
                collision += count;
            }
        }

        let mut hist = histogram::Config::new().precision(5).build().unwrap();
        for i in 0..(index.sa.offsets.len() - 1) {
            let len = index.sa.offsets[i + 1] - index.sa.offsets[i];
            hist.increment(len as u64).unwrap();
        }

        // k buckets kmers util_count coll_count util coll max p50 p90 p99 p999 stddev
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            index.sa.bucket_width,
            index.sa.offsets.len() - 1,
            kmers.len(),
            non_zero,
            collision,
            non_zero as f64 * 100.0 / (index.sa.offsets.len() - 1) as f64,
            collision as f64 * 100.0 / kmers.len() as f64,
            hist.maximum().unwrap(),
            hist.percentile(50.0).unwrap(),
            hist.percentile(90.0).unwrap(),
            hist.percentile(99.0).unwrap(),
            hist.percentile(99.9).unwrap(),
            hist.stddev().unwrap(),
        );

        Ok(())
    }
}
