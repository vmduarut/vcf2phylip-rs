use crate::matrix::Matrix;
use anyhow::{Result, ensure};
use std::io::{BufWriter, Write};

/// Sequential relaxed PHYLIP: full identifiers, two-space separator, no truncation.
pub fn write(writer: impl Write, matrix: &Matrix, order: &[usize]) -> Result<()> {
    ensure!(matrix.sites() > 0, "cannot write a zero-site alignment");
    matrix.validate_order(order)?;
    let mut writer = BufWriter::new(writer);
    writeln!(writer, "{} {}", matrix.samples().len(), matrix.sites())?;
    for &i in order {
        write!(writer, "{}  ", matrix.samples()[i])?;
        writer.write_all(&matrix.sequences()[i])?;
        writer.write_all(b"\n")?;
    }
    writer.flush()?;
    Ok(())
}
