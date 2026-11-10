use crate::matrix::Alignment;
use anyhow::{Result, ensure};
use std::io::{BufWriter, Write};

/// Sequential relaxed PHYLIP: full identifiers, two-space separator, no truncation.
pub fn write(
    writer: impl Write,
    matrix: &(impl Alignment + ?Sized),
    order: &[usize],
) -> Result<()> {
    ensure!(matrix.sites() > 0, "cannot write a zero-site alignment");
    matrix.validate_order(order)?;
    let mut writer = BufWriter::new(writer);
    writeln!(writer, "{} {}", matrix.samples().len(), matrix.sites())?;
    for &i in order {
        write!(writer, "{}  ", matrix.samples()[i])?;
        matrix.row_chunks(i, &mut |chunk| Ok(writer.write_all(chunk)?))?;
        writer.write_all(b"\n")?;
    }
    writer.flush()?;
    Ok(())
}
