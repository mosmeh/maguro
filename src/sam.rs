use crate::index::{Index, SequenceId};
use std::io;

pub struct SamWriter<W: io::Write> {
    out: W,
}

impl SamWriter<io::Stdout> {
    pub fn from_stdout() -> SamWriter<io::Stdout> {
        SamWriter { out: io::stdout() }
    }
}

impl<W: io::Write> SamWriter<W> {
    pub fn write_header(&mut self, index: &Index) -> io::Result<()> {
        self.out.write_all(b"@HD\tVN:1.0\tSO:unknown\n")?;

        for i in 0..index.num_seqs() {
            let id = SequenceId(i);
            self.out.write_all(b"@SQ\tSN:")?;
            self.out.write_all(index.seq_name(id))?;
            let range = index.seq_range(id);
            writeln!(self.out, "\tLN:{}", range.end - range.start)?;
        }

        writeln!(
            self.out,
            "@PG\tID:{}\tPN:{}\tVN:{}",
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_VERSION")
        )
    }

    pub fn write_mapping(
        &mut self,
        qname: &[u8],
        rname: &[u8],
        seq: &[u8],
        pos: usize,
        rc: bool,
    ) -> io::Result<()> {
        let flag = if rc { 0x10 } else { 0 };

        self.out.write_all(qname)?;
        write!(self.out, "\t{}\t", flag)?;
        self.out.write_all(rname)?;
        write!(self.out, "\t{}\t255\t{}M\t*\t0\t0\t", pos + 1, seq.len())?;
        self.out.write_all(seq)?;
        self.out.write_all(b"\t*\n")
    }

    pub fn write_unmapped(&mut self, qname: &[u8], seq: &[u8]) -> io::Result<()> {
        self.out.write_all(qname)?;
        self.out.write_all(b"\t4\t*\t0\t255\t*\t*\t0\t0\t")?;
        self.out.write_all(seq)?;
        self.out.write_all(b"\t*\n")
    }
}
