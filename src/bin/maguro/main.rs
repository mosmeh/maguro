mod commands;

use commands::*;
use structopt::{clap::AppSettings, StructOpt};

#[cfg(all(target_env = "musl", target_pointer_width = "64"))]
#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

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
    Stats(StatsCommand),
}

fn main() -> anyhow::Result<()> {
    match Opt::from_args() {
        Opt::Index(cmd) => cmd.run(),
        Opt::Map(cmd) => cmd.run(),
        Opt::Stats(cmd) => cmd.run(),
    }
}
