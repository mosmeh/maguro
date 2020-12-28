use crate::{
    index::{Index, SequenceId},
    mapper::Mapping,
};
use std::io;

pub struct SamWriter<W: io::Write> {
    out: W,
}

impl<W: io::Write> SamWriter<W> {
    pub fn new(out: W) -> Self {
        SamWriter { out }
    }

    #[allow(clippy::write_with_newline)]
    pub fn write_header(&mut self, index: &Index) -> io::Result<()> {
        self.out.write_all(b"@HD\tVN:1.0\tSO:unknown\n")?;

        for i in 0..index.num_seqs() {
            let id = SequenceId(i);
            self.out.write_all(b"@SQ\tSN:")?;
            self.out.write_all(index.seq_name(id))?;
            let range = index.seq_range(id);
            write!(self.out, "\tLN:{}\n", range.end - range.start)?;
        }

        write!(
            self.out,
            "@PG\tID:{}\tPN:{}\tVN:{}\n",
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_VERSION")
        )
    }

    pub fn write_mapping(
        &mut self,
        index: &Index,
        qname: &[u8],
        seq: &[u8],
        mapping: &Mapping,
        secondary: bool,
    ) -> io::Result<()> {
        let flag = if secondary { 0x100 } else { 0 };
        let rname = index.seq_name(mapping.seq_id);

        self.out.write_all(qname)?;
        write!(self.out, "\t{}\t", flag)?;
        self.out.write_all(rname)?;
        write!(
            self.out,
            "\t{}\t255\t{}M\t*\t0\t{}\t",
            mapping.pos + 1,
            seq.len(),
            seq.len()
        )?;
        self.out.write_all(&seq)?;
        self.out.write_all(b"\t*\n")
    }

    pub fn write_unmapped(&mut self, qname: &[u8], seq: &[u8]) -> io::Result<()> {
        self.out.write_all(qname)?;
        self.out.write_all(b"\t4\t*\t0\t255\t*\t*\t0\t0\t")?;
        self.out.write_all(seq)?;
        self.out.write_all(b"\t*\n")
    }
}
