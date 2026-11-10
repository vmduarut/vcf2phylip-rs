//! No-clobber staging. Finalization is atomic per file, not per output set.
use anyhow::{Context, Result, bail, ensure};
use std::{
    fs,
    path::{Path, PathBuf},
};
use tempfile::NamedTempFile;

pub(crate) fn prefix(input: &Path, explicit: Option<&str>) -> Result<String> {
    let name = if let Some(name) = explicit {
        name.to_owned()
    } else {
        let mut name = input
            .file_name()
            .and_then(|s| s.to_str())
            .context("input basename is not UTF-8; supply --output-prefix")?
            .to_owned();
        let lower = name.to_ascii_lowercase();
        for suffix in [
            ".vcf.gz",
            ".vcf.bgz",
            ".vcf.bgzf",
            ".bcf.gz",
            ".bcf.bgz",
            ".bcf.bgzf",
            ".vcf",
            ".bcf",
        ] {
            if lower.ends_with(suffix) {
                name.truncate(name.len() - suffix.len());
                break;
            }
        }
        name
    };
    ensure!(
        !name.is_empty()
            && name != "."
            && name != ".."
            && !name.contains(['/', '\\'])
            && !name.chars().any(char::is_control),
        "output prefix must be a nonempty basename, not a path"
    );
    Ok(name)
}

pub(crate) fn preflight(paths: &[PathBuf]) -> Result<()> {
    for path in paths {
        match fs::symlink_metadata(path) {
            Ok(_) => bail!(
                "refusing existing output target (no overwrite): {}",
                path.display()
            ),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e).with_context(|| format!("checking output {}", path.display())),
        }
    }
    Ok(())
}

pub(crate) struct Staged {
    pub file: NamedTempFile,
    pub path: PathBuf,
}

impl Staged {
    pub fn new(path: PathBuf) -> Result<Self> {
        let file = NamedTempFile::new_in(path.parent().context("output has no parent")?)?;
        Ok(Self { file, path })
    }

    pub fn persist(self) -> Result<PathBuf> {
        self.file.persist_noclobber(&self.path).with_context(|| {
            format!(
                "finalizing {}; earlier outputs may already have been finalized",
                self.path.display()
            )
        })?;
        Ok(self.path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn racing_creator_is_never_overwritten_and_staging_is_cleaned() {
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("alignment.phy");
        preflight(std::slice::from_ref(&target)).unwrap();
        let mut staged = Staged::new(target.clone()).unwrap();
        staged.file.write_all(b"our alignment").unwrap();
        fs::write(&target, b"racing creator").unwrap();
        assert!(staged.persist().is_err());
        assert_eq!(fs::read(&target).unwrap(), b"racing creator");
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    #[test]
    fn prefix_suffixes_and_paths() {
        for (input, expected) in [
            ("dir/name.VCF.GZ", "name"),
            ("name.BCF", "name"),
            ("name.vcf.middle.bcf.bgzf", "name.vcf.middle"),
            ("name.other.gz", "name.other.gz"),
        ] {
            assert_eq!(prefix(Path::new(input), None).unwrap(), expected);
        }
        for invalid in ["", ".", "..", "../x", "/x", "a/b", "a\\b", "a\n"] {
            assert!(prefix(Path::new("input"), Some(invalid)).is_err());
        }
        assert_eq!(
            prefix(Path::new("input"), Some("study..v1")).unwrap(),
            "study..v1"
        );
    }
}
