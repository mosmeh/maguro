mod index;
mod map;
mod stats;

pub use index::IndexCommand;
pub use map::MapCommand;
pub use stats::StatsCommand;

pub trait Command {
    fn run(self) -> anyhow::Result<()>;
}
