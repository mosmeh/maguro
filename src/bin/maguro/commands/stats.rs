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

        /*let mut non_zero = 0;
        for i in 0..(index.sa.offsets.len() - 1) {
            let len = index.sa.offsets[i + 1] - index.sa.offsets[i];
            if len > 0 {
                non_zero += 1;
            }
        }

        let mut hist = histogram::Config::new().precision(5).build().unwrap();
        for i in 0..(index.sa.offsets.len() - 1) {
            let len = index.sa.offsets[i + 1] - index.sa.offsets[i];
            hist.increment(len as u64).unwrap();
        }*/

        // k f array offsets buckets
        println!(
            "{}\t{}\t{}\t{}\t{}",
            index.sa.k,
            index.sa.f,
            index.sa.array.len(),
            index.sa.offsets.len(),
            index.sa.buckets.len()
        );

        Ok(())
    }
}
