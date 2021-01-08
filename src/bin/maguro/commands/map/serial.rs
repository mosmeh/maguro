use super::MapCommand;
use anyhow::{anyhow, Result};
use bio::io::fasta::{self, FastaRead};
use maguro::{index::Index, mapper::Mapper, utils};
use std::io::{self, BufWriter, Write};

pub fn main(config: MapCommand, index: &Index, mapper: &Mapper) -> Result<()> {
    let out = io::stdout();
    let mut out = BufWriter::new(out.lock());

    let mut num_processed = 0;
    let mut num_mapped = 0;

    if let Some(ref reads) = config.reads {
        eprintln!("Starting single-end mapping");

        let mut reader = fasta::Reader::from_file(&reads)?;
        let mut record = fasta::Record::new();
        reader.read(&mut record)?;

        while !record.is_empty() {
            record.check().map_err(|e| anyhow!(e.to_owned()))?;

            let qname = utils::extract_name_bytes(record.id(), &config.header_sep);
            let mapped = super::map_single(&mut out, &index, &mapper, qname, record.seq())?;

            if mapped {
                num_mapped += 1
            };
            num_processed += 1;

            reader.read(&mut record)?;
        }
    } else {
        eprintln!("Starting paired-end mapping");

        let mut reader1 = fasta::Reader::from_file(config.mates1.as_ref().unwrap())?;
        let mut reader2 = fasta::Reader::from_file(config.mates2.as_ref().unwrap())?;

        let mut record1 = fasta::Record::new();
        let mut record2 = fasta::Record::new();
        reader1.read(&mut record1)?;
        reader2.read(&mut record2)?;

        while !record1.is_empty() {
            record1.check().map_err(|e| anyhow!(e.to_owned()))?;
            record2.check().map_err(|e| anyhow!(e.to_owned()))?;

            let qname = utils::extract_name_bytes(record1.id(), &config.header_sep);
            let mapped = super::map_pair(
                &mut out,
                &index,
                &mapper,
                qname,
                record1.seq(),
                record2.seq(),
            )?;

            if mapped {
                num_mapped += 1
            };
            num_processed += 1;

            reader1.read(&mut record1)?;
            reader2.read(&mut record2)?;
        }

        assert!(record2.is_empty());
    }

    out.flush()?;

    eprintln!(
        "Mapped {} / {} reads ({:.2}%)",
        num_mapped,
        num_processed,
        num_mapped as f64 * 100.0 / num_processed as f64
    );

    Ok(())
}
