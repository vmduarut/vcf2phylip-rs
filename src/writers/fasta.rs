use crate::matrix::Matrix;
use anyhow::{Result, ensure};
use std::io::{BufWriter, Write};

/// Identifiers are preserved verbatim; sequence lines wrap at 80 bytes.
pub fn write(writer: impl Write, matrix: &Matrix, order: &[usize]) -> Result<()> {
    ensure!(matrix.sites() > 0, "cannot write a zero-site alignment");
    matrix.validate_order(order)?;
    let mut writer = BufWriter::new(writer);
    for &i in order {
        writeln!(writer, ">{}", matrix.samples()[i])?;
        for chunk in matrix.sequences()[i].chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}
