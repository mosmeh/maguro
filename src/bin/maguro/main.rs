mod commands;

use commands::{Command, IndexCommand, MapCommand};
use structopt::{clap::AppSettings, StructOpt};

#[derive(StructOpt)]
#[structopt(
    name = env!("CARGO_PKG_NAME"),
    author = env!("CARGO_PKG_AUTHORS"),
    about = env!("CARGO_PKG_DESCRIPTION"),
    rename_all = "kebab-case",
    setting(AppSettings::ColoredHelp),
    setting(AppSettings::DeriveDisplayOrder),
    setting(AppSettings::AllArgsOverrideSelf)
)]
enum Opt {
    Index(IndexCommand),
    Map(MapCommand),
}

fn main() -> anyhow::Result<()> {
    match Opt::from_args() {
        Opt::Index(cmd) => cmd.run(),
        Opt::Map(cmd) => cmd.run(),
    }
}
