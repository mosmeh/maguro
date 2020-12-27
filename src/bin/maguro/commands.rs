mod index;
mod map;

pub use index::IndexCommand;
pub use map::MapCommand;

pub trait Command {
    fn run(self) -> anyhow::Result<()>;
}
