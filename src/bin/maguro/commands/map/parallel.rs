use super::MapCommand;
use anyhow::{anyhow, Result};
use bio::io::fasta::{self, FastaRead};
use crossbeam_deque::{Injector, Stealer, Worker};
use maguro::{index::Index, mapper::Mapper, utils};
use std::io::{Stdout, Write};
use std::sync::atomic::{AtomicBool, Ordering};

pub fn parallel_map(
    config: MapCommand,
    out: &mut Stdout,
    index: &Index,
    reader: &mut fasta::Reader<std::fs::File>,
) -> Result<()> {
    let num_workers = config.threads.saturating_sub(1).max(1);
    let input_finished = AtomicBool::new(false);
    let injector = Injector::new();
    let workers: Vec<_> = (0..num_workers).map(|_| Worker::new_fifo()).collect();
    let stealers: Vec<_> = workers.iter().map(|worker| worker.stealer()).collect();

    crossbeam_utils::thread::scope(|s| -> Result<()> {
        for (worker_id, worker) in workers.into_iter().enumerate() {
            let stealers: Vec<_> = stealers
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != worker_id)
                .map(|(_, stealer)| stealer.clone())
                .collect();
            s.spawn(|_| {
                MapWorker::new(
                    &config,
                    &out,
                    &index,
                    &input_finished,
                    &injector,
                    worker,
                    stealers,
                )
                .run()
                .unwrap()
            });
        }

        let mut record = fasta::Record::new();
        reader.read(&mut record)?;

        while !record.is_empty() {
            record.check().map_err(|e| anyhow!(e.to_owned()))?;

            injector.push(Task {
                qname: utils::extract_name_bytes(record.id(), &config.header_sep).to_owned(),
                seq: record.seq().to_owned(),
            });

            reader.read(&mut record)?;
        }

        input_finished.store(true, Ordering::Relaxed);

        Ok(())
    })
    .unwrap()
}

struct Task {
    qname: Vec<u8>,
    seq: Vec<u8>,
}

struct MapWorker<'a> {
    out: &'a Stdout,
    index: &'a Index,
    input_finished: &'a AtomicBool,
    injector: &'a Injector<Task>,
    worker: Worker<Task>,
    stealers: Vec<Stealer<Task>>,
    mapper: Mapper<'a>,
    out_buf: Vec<u8>,
}

impl<'a> MapWorker<'a> {
    fn new(
        config: &MapCommand,
        out: &'a Stdout,
        index: &'a Index,
        input_finished: &'a AtomicBool,
        injector: &'a Injector<Task>,
        worker: Worker<Task>,
        stealers: Vec<Stealer<Task>>,
    ) -> Self {
        Self {
            out,
            index,
            mapper: super::build_mapper(index, config),
            injector,
            worker,
            stealers,
            input_finished,
            out_buf: Vec::new(),
        }
    }

    fn find_task(&mut self) -> Option<Task> {
        self.worker.pop().or_else(|| {
            std::iter::repeat_with(|| {
                self.injector.steal_batch_and_pop(&self.worker).or_else(|| {
                    self.stealers
                        .iter()
                        .map(|stealer| stealer.steal())
                        .collect()
                })
            })
            .find(|steal| !steal.is_retry())
            .and_then(|steal| steal.success())
        })
    }

    fn run(&mut self) -> Result<()> {
        loop {
            while let Some(task) = self.find_task() {
                super::map(
                    &mut self.out_buf,
                    self.index,
                    &mut self.mapper,
                    &task.qname,
                    &task.seq,
                )?;
            }

            if !self.out_buf.is_empty() {
                {
                    let mut out = self.out.lock();
                    out.write_all(&self.out_buf)?;
                }
                self.out_buf.clear();
            }

            if self.input_finished.load(Ordering::Relaxed) && self.injector.is_empty() {
                return Ok(());
            }
        }
    }
}
