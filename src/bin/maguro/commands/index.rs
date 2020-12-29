use super::Command;
use maguro::index::IndexBuilder;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};
use structopt::StructOpt;

#[derive(StructOpt)]
pub struct IndexCommand {
    #[structopt(short, long)]
    reference: PathBuf,
    #[structopt(short, long)]
    output: PathBuf,
    #[structopt(short, long, default_value = "10")]
    bucket_width: usize,
    #[structopt(long)]
    header_sep: Option<String>,
}

impl Command for IndexCommand {
    fn run(self) -> anyhow::Result<()> {
        let mut builder = IndexBuilder::from_file(self.reference)?.bucket_width(self.bucket_width);
        if let Some(value) = self.header_sep {
            builder = builder.header_sep(value);
        }

        let index = builder.build()?;

        let mut writer = BufWriter::new(File::create(&self.output)?);
        bincode::serialize_into(&mut writer, &index)?;
        writer.flush()?;

        Ok(())
    }
}
