use super::Command;
use maguro::index::IndexBuilder;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
pub struct IndexCommand {
    #[structopt(short, long)]
    reference: PathBuf,
    #[structopt(short, long)]
    index: PathBuf,
    #[structopt(short, default_value = "13")]
    k: usize,
    #[structopt(short, default_value = "1")]
    f: f64,
    #[structopt(long)]
    header_sep: Option<String>,
}

impl Command for IndexCommand {
    fn run(self) -> anyhow::Result<()> {
        eprintln!("{:#?}", self);

        let mut builder = IndexBuilder::from_file(self.reference)?.k(self.k).f(self.f);
        if let Some(value) = self.header_sep {
            builder = builder.header_sep(value);
        }

        eprintln!("Indexing");
        let index = builder.build()?;

        eprintln!("Writing");
        let mut writer = BufWriter::new(File::create(&self.index)?);
        bincode::serialize_into(&mut writer, &index)?;
        writer.flush()?;

        Ok(())
    }
}
