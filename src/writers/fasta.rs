use crate::matrix::Alignment;
use anyhow::{Result, ensure};
use std::io::{BufWriter, Write};

/// Identifiers are preserved verbatim; sequence lines wrap at 80 bytes.
pub fn write(
    writer: impl Write,
    matrix: &(impl Alignment + ?Sized),
    sample_order: &[usize],
) -> Result<()> {
    ensure!(matrix.sites() > 0, "cannot write a zero-site alignment");
    matrix.validate_order(sample_order)?;
    let mut writer = BufWriter::new(writer);
    for &i in sample_order {
        writeln!(writer, ">{}", matrix.samples()[i])?;
        let mut column = 0;
        matrix.row_chunks(i, &mut |mut chunk| {
            while !chunk.is_empty() {
                let n = (80 - column).min(chunk.len());
                writer.write_all(&chunk[..n])?;
                column += n;
                chunk = &chunk[n..];
                if column == 80 {
                    writer.write_all(b"\n")?;
                    column = 0;
                }
            }
            Ok(())
        })?;
        if column != 0 {
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}
