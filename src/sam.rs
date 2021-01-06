#![allow(clippy::write_with_newline)]

use crate::{
    index::{Index, SequenceId},
    mapper::{pair::PairMapping, single::SingleMapping},
    sequence,
};
use std::io;

pub fn write_header<W: io::Write>(mut out: W, index: &Index) -> io::Result<()> {
    out.write_all(b"@HD\tVN:1.0\tSO:unknown\n")?;

    for i in 0..index.num_seqs() {
        let id = SequenceId(i);
        out.write_all(b"@SQ\tSN:")?;
        out.write_all(index.seq_name(id))?;
        let range = index.seq_range(id);
        write!(out, "\tLN:{}\n", range.end - range.start)?;
    }

    write!(
        out,
        "@PG\tID:{}\tPN:{}\tVN:{}\n",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION")
    )
}

pub fn write_mapping_single<W: io::Write>(
    mut out: W,
    index: &Index,
    qname: &[u8],
    seq: &[u8],
    mapping: &SingleMapping,
    secondary: bool,
    rc_seq_cache: &mut Option<Vec<u8>>,
) -> io::Result<()> {
    let rname = index.seq_name(mapping.seq_id);

    let mut flag = 0;
    if mapping.strand.is_reverse() {
        flag |= 0x10;
    }
    if secondary {
        flag |= 0x100;
    }

    let seq: &[u8] = if mapping.strand.is_forward() {
        seq
    } else if let Some(rc) = rc_seq_cache {
        rc
    } else {
        *rc_seq_cache = Some(sequence::reverse_complement(seq));
        &rc_seq_cache.as_ref().unwrap()
    };

    out.write_all(qname)?;
    write!(out, "\t{}\t", flag)?;
    out.write_all(rname)?;
    write!(out, "\t{}\t255\t{}M\t*\t0\t0\t", mapping.pos + 1, seq.len())?;
    out.write_all(seq)?;
    write!(out, "\t*\tAS:i:{}\n", mapping.score)
}

pub fn write_unmapped_single<W: io::Write>(mut out: W, qname: &[u8], seq: &[u8]) -> io::Result<()> {
    out.write_all(qname)?;
    out.write_all(b"\t4\t*\t0\t255\t*\t*\t0\t0\t")?;
    out.write_all(seq)?;
    out.write_all(b"\t*\n")
}

#[allow(clippy::too_many_arguments)]
pub fn write_mapping_pair<W: io::Write>(
    mut out: W,
    index: &Index,
    qname: &[u8],
    seq1: &[u8],
    seq2: &[u8],
    mapping: &PairMapping,
    secondary: bool,
    rc_seq_cache1: &mut Option<Vec<u8>>,
    rc_seq_cache2: &mut Option<Vec<u8>>,
) -> io::Result<()> {
    let rname = index.seq_name(mapping.seq_id);

    let mut flag1 = 0x1 | 0x2 | 0x40;
    let mut flag2 = 0x1 | 0x2 | 0x80;
    if mapping.strand.is_forward() {
        flag1 |= 0x20;
        flag2 |= 0x10;
    } else {
        flag1 |= 0x10;
        flag2 |= 0x20;
    }
    if secondary {
        flag1 |= 0x100;
        flag2 |= 0x100;
    }

    let seq1: &[u8] = if mapping.strand.is_forward() {
        seq1
    } else if let Some(rc) = rc_seq_cache1 {
        rc
    } else {
        *rc_seq_cache1 = Some(sequence::reverse_complement(seq1));
        &rc_seq_cache1.as_ref().unwrap()
    };
    let seq2: &[u8] = if mapping.strand.is_reverse() {
        seq2
    } else if let Some(rc) = rc_seq_cache2 {
        rc
    } else {
        *rc_seq_cache2 = Some(sequence::reverse_complement(seq2));
        &rc_seq_cache2.as_ref().unwrap()
    };

    // TODO: handle qname1 != qname2

    out.write_all(qname)?;
    write!(out, "\t{}\t", flag1)?;
    out.write_all(rname)?;
    write!(
        out,
        "\t{}\t255\t{}M\t*\t0\t0\t",
        mapping.pos1 + 1,
        seq1.len()
    )?;
    out.write_all(seq1)?;
    out.write_all(b"\t*\n")?;

    out.write_all(qname)?;
    write!(out, "\t{}\t", flag2)?;
    out.write_all(rname)?;
    write!(
        out,
        "\t{}\t255\t{}M\t*\t0\t0\t",
        mapping.pos2 + 1,
        seq2.len()
    )?;
    out.write_all(seq2)?;
    write!(out, "\t*\tAS:i:{}\n", mapping.score)
}

pub fn write_unmapped_pair<W: io::Write>(
    mut out: W,
    qname: &[u8],
    seq1: &[u8],
    seq2: &[u8],
) -> io::Result<()> {
    let flag = 0x1 | 0x4 | 0x8 | 0x40;
    out.write_all(qname)?;
    write!(out, "\t{}\t*\t0\t255\t*\t*\t0\t0\t", flag)?;
    out.write_all(seq1)?;
    out.write_all(b"\t*\n")?;

    let flag = 0x1 | 0x4 | 0x8 | 0x80;
    out.write_all(qname)?;
    write!(out, "\t{}\t*\t0\t255\t*\t*\t0\t0\t", flag)?;
    out.write_all(seq2)?;
    out.write_all(b"\t*\n")
}
