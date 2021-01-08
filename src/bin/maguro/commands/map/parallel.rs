use super::MapCommand;
use anyhow::{anyhow, Result};
use bio::io::fasta::{self, FastaRead};
use maguro::{index::Index, mapper::Mapper, utils};
use rayon::prelude::*;
use std::{
    io::{self, BufWriter, Write},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
    thread,
};

struct SingleTask {
    qname: Vec<u8>,
    seq: Vec<u8>,
}

struct PairTask {
    qname: Vec<u8>,
    seq1: Vec<u8>,
    seq2: Vec<u8>,
}

pub fn main(config: MapCommand, index: &Index, mapper: &Mapper) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()?;

    let chunk_size = config.chunk * 1024 * 1024;

    let (writer_tx, writer_rx): (crossbeam_channel::Sender<Vec<u8>>, _) =
        crossbeam_channel::unbounded();
    let writer_thread = thread::spawn(move || -> Result<()> {
        let out = io::stdout();
        let mut writer = BufWriter::new(out.lock());
        for x in writer_rx.into_iter() {
            writer.write_all(&x)?;
        }
        Ok(())
    });

    let mut num_processed = 0;
    let num_mapped = Arc::new(AtomicUsize::new(0));

    match (config.reads, config.mates1, config.mates2) {
        (Some(reads), None, None) => {
            eprintln!("Starting single-end mapping");

            let mut reader = fasta::Reader::from_file(reads)?;
            let mut record = fasta::Record::new();
            reader.read(&mut record)?;

            while !record.is_empty() {
                let mut chunk = Vec::new();

                while !record.is_empty() && chunk.len() < chunk_size {
                    record.check().map_err(|e| anyhow!(e.to_owned()))?;
                    chunk.push(SingleTask {
                        qname: utils::extract_name_bytes(record.id(), &config.header_sep)
                            .to_owned(),
                        seq: record.seq().to_owned(),
                    });
                    reader.read(&mut record)?;
                }

                chunk.par_iter().try_for_each_with::<_, _, Result<()>>(
                    writer_tx.clone(),
                    |tx, task| {
                        let mut buf = Vec::new();
                        let mapped =
                            super::map_single(&mut buf, index, mapper, &task.qname, &task.seq)?;
                        if mapped {
                            num_mapped.fetch_add(1, Ordering::Relaxed);
                        }
                        tx.send(buf)?;
                        Ok(())
                    },
                )?;

                num_processed += chunk.len();
            }
        }
        (None, Some(mates1), Some(mates2)) => {
            eprintln!("Starting paired-end mapping");

            let mut reader1 = fasta::Reader::from_file(mates1)?;
            let mut reader2 = fasta::Reader::from_file(mates2)?;

            let mut record1 = fasta::Record::new();
            let mut record2 = fasta::Record::new();
            reader1.read(&mut record1)?;
            reader2.read(&mut record2)?;

            while !record1.is_empty() {
                let mut chunk = Vec::new();

                while !record1.is_empty() && chunk.len() < chunk_size {
                    record1.check().map_err(|e| anyhow!(e.to_owned()))?;
                    record2.check().map_err(|e| anyhow!(e.to_owned()))?;
                    chunk.push(PairTask {
                        qname: utils::extract_name_bytes(record1.id(), &config.header_sep)
                            .to_owned(),
                        seq1: record1.seq().to_owned(),
                        seq2: record2.seq().to_owned(),
                    });
                    reader1.read(&mut record1)?;
                    reader2.read(&mut record2)?;
                }

                chunk.par_iter().try_for_each_with::<_, _, Result<()>>(
                    (writer_tx.clone(), num_mapped.clone()),
                    |(tx, num_mapped), task| {
                        let mut buf = Vec::new();
                        let mapped = super::map_pair(
                            &mut buf,
                            index,
                            mapper,
                            &task.qname,
                            &task.seq1,
                            &task.seq2,
                        )?;
                        if mapped {
                            num_mapped.fetch_add(1, Ordering::Relaxed);
                        }
                        tx.send(buf)?;
                        Ok(())
                    },
                )?;

                num_processed += chunk.len();
            }

            assert!(record2.is_empty());
        }
        _ => unreachable!(),
    }

    eprintln!("Finishing output");
    drop(writer_tx);
    writer_thread.join().unwrap()?;

    let num_mapped = num_mapped.load(std::sync::atomic::Ordering::Relaxed);
    eprintln!(
        "Mapped {} / {} reads ({:.2}%)",
        num_mapped,
        num_processed,
        num_mapped as f64 * 100.0 / num_processed as f64
    );

    Ok(())
}
