//! Uniform sample-major tiles in one managed spool; no per-locus/tile index.
use crate::{
    filter::ProcessedLocus,
    matrix::{Alignment, Matrix},
};
use anyhow::{Context, Result, ensure};
use std::{
    cell::RefCell,
    io::{Read, Seek, SeekFrom, Write},
    path::Path,
};
use tempfile::NamedTempFile;

#[derive(Clone, Copy, Debug, Default, clap::ValueEnum, PartialEq, Eq)]
pub enum StorageMode {
    #[default]
    Memory,
    Disk,
}

/// Bounded-cell-buffer alignment builder. Call `finish` before reading.
/// Any append/finish I/O error is terminal: discard the builder rather than retry.
pub struct DiskMatrix {
    samples: Vec<String>,
    file: RefCell<NamedTempFile>,
    buffer: RefCell<Vec<u8>>,
    width: usize,
    pending: usize,
    sites: usize,
    finished: bool,
}

impl DiskMatrix {
    /// Budget covers the cell buffer's actual capacity only (not process RSS).
    /// At least one byte per sample is required. Header/record state is extra.
    pub fn new(samples: Vec<String>, budget: usize, temp_dir: Option<&Path>) -> Result<Self> {
        Matrix::new(samples.clone())?;
        ensure!(
            budget >= samples.len(),
            "--memory-limit must be at least the sample count ({} bytes)",
            samples.len()
        );
        // Avoid allocating needlessly large tiles for small sample counts.
        let width = (budget / samples.len()).min(65536);
        let size = width
            .checked_mul(samples.len())
            .context("tile size overflow")?;
        let mut buffer = Vec::new();
        buffer.try_reserve_exact(size)?;
        ensure!(
            buffer.capacity() <= budget,
            "allocator tile capacity exceeds --memory-limit"
        );
        buffer.resize(size, 0);
        let file = match temp_dir {
            Some(dir) => NamedTempFile::new_in(dir),
            None => NamedTempFile::new(),
        }
        .context("creating matrix spool in temporary directory")?;
        Ok(Self {
            samples,
            file: RefCell::new(file),
            buffer: RefCell::new(buffer),
            width,
            pending: 0,
            sites: 0,
            finished: false,
        })
    }

    pub fn append(&mut self, locus: &ProcessedLocus) -> Result<()> {
        ensure!(!self.finished, "cannot append to finished spool");
        if !locus.nucleotide_accepted {
            return Ok(());
        }
        ensure!(
            locus.nucleotide.len() == self.samples.len(),
            "matrix/sample dimension mismatch"
        );
        let sites = self
            .sites
            .checked_add(1)
            .context("matrix length overflow")?;
        let bytes = sites
            .checked_mul(self.samples.len())
            .context("matrix size overflow")?;
        u64::try_from(bytes).context("matrix file offset overflow")?;
        for (sample, call) in locus.nucleotide.iter().enumerate() {
            self.buffer.get_mut()[sample * self.width + self.pending] = call.symbol();
        }
        self.pending += 1;
        self.sites = sites;
        if self.pending == self.width {
            self.flush_tile()?;
        }
        Ok(())
    }

    fn flush_tile(&mut self) -> Result<()> {
        if self.pending == 0 {
            return Ok(());
        }
        for row in self.buffer.get_mut().chunks(self.width) {
            self.file
                .get_mut()
                .write_all(&row[..self.pending])
                .context("writing matrix spool")?;
        }
        self.pending = 0;
        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        self.flush_tile()?;
        self.file.get_mut().flush()?;
        self.finished = true;
        Ok(())
    }
}

impl Alignment for DiskMatrix {
    fn samples(&self) -> &[String] {
        &self.samples
    }
    fn sites(&self) -> usize {
        self.sites
    }
    fn row_chunks(&self, sample: usize, emit: &mut dyn FnMut(&[u8]) -> Result<()>) -> Result<()> {
        ensure!(self.finished, "finish spool before reading");
        ensure!(sample < self.samples.len(), "sample index out of bounds");
        let mut file = self.file.try_borrow_mut()?;
        let mut buffer = self.buffer.try_borrow_mut()?;
        for start in (0..self.sites).step_by(self.width) {
            let width = self.width.min(self.sites - start);
            let offset = start
                .checked_mul(self.samples.len())
                .and_then(|x| sample.checked_mul(width).and_then(|y| x.checked_add(y)))
                .context("matrix offset overflow")?;
            file.seek(SeekFrom::Start(u64::try_from(offset)?))?;
            file.read_exact(&mut buffer[..width])
                .context("reading matrix spool")?;
            emit(&buffer[..width])?;
        }
        Ok(())
    }
}

pub enum Storage {
    Memory(Matrix),
    Disk(DiskMatrix),
}
impl Storage {
    pub fn new(
        samples: Vec<String>,
        mode: StorageMode,
        budget: usize,
        temp: Option<&Path>,
    ) -> Result<Self> {
        Ok(match mode {
            StorageMode::Memory => Self::Memory(Matrix::new(samples)?),
            StorageMode::Disk => Self::Disk(DiskMatrix::new(samples, budget, temp)?),
        })
    }
    pub fn append(&mut self, locus: &ProcessedLocus) -> Result<()> {
        match self {
            Self::Memory(m) => m.append(locus),
            Self::Disk(m) => m.append(locus),
        }
    }
    pub fn finish(&mut self) -> Result<()> {
        match self {
            Self::Memory(_) => Ok(()),
            Self::Disk(m) => m.finish(),
        }
    }
}
impl Alignment for Storage {
    fn samples(&self) -> &[String] {
        match self {
            Self::Memory(m) => m.samples(),
            Self::Disk(m) => m.samples(),
        }
    }
    fn sites(&self) -> usize {
        match self {
            Self::Memory(m) => m.sites(),
            Self::Disk(m) => m.sites(),
        }
    }
    fn row_chunks(&self, sample: usize, emit: &mut dyn FnMut(&[u8]) -> Result<()>) -> Result<()> {
        match self {
            Self::Memory(m) => m.row_chunks(sample, emit),
            Self::Disk(m) => m.row_chunks(sample, emit),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        filter::ProcessingConfig,
        genotype::{PartialMissingPolicy, Resolution},
        reader::VariantReader,
    };
    use std::io::Cursor;

    fn locus() -> ProcessedLocus {
        let mut reader = VariantReader::from_reader(Cursor::new(b"##fileformat=VCFv4.3\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n")).unwrap();
        ProcessingConfig::new(1, Some(1), PartialMissingPolicy::Missing, Resolution::Iupac)
            .unwrap()
            .process(&reader.next_variant().unwrap().unwrap(), 0)
            .unwrap()
    }

    #[test]
    fn disk_read_callback_and_prewrite_errors() {
        let dir = tempfile::tempdir().unwrap();
        let mut matrix = DiskMatrix::new(vec!["S".into()], 1, Some(dir.path())).unwrap();
        let mut output = Vec::new();
        assert!(crate::writers::fasta::write(&mut output, &matrix, &[0]).is_err());
        assert!(crate::writers::phylip::write(&mut output, &matrix, &[0]).is_err());
        assert!(output.is_empty());
        matrix.append(&locus()).unwrap();
        matrix.finish().unwrap();
        assert!(crate::writers::fasta::write(&mut output, &matrix, &[1]).is_err());
        assert!(output.is_empty());
        assert!(
            matrix
                .row_chunks(0, &mut |_| anyhow::bail!("injected callback failure"))
                .is_err()
        );
        matrix.file.get_mut().as_file_mut().set_len(0).unwrap();
        assert!(matrix.row_chunks(0, &mut |_| Ok(())).is_err());
        drop(matrix);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn disk_full_write_propagates_and_cleans_spool() {
        let dir = tempfile::tempdir().unwrap();
        let mut matrix = DiskMatrix::new(vec!["S".into()], 1, Some(dir.path())).unwrap();
        // Replace only the open handle; NamedTempFile still owns its original path.
        *matrix.file.get_mut().as_file_mut() = std::fs::OpenOptions::new()
            .write(true)
            .open("/dev/full")
            .unwrap();
        assert!(
            matrix
                .append(&locus())
                .unwrap_err()
                .to_string()
                .contains("writing matrix spool")
        );
        drop(matrix);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }
}
