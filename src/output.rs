//! Output staging. Finalization is atomic per file, not per output set.
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

pub(crate) fn preflight(paths: &[PathBuf], input: &Path, force: bool) -> Result<()> {
    for path in paths {
        check_target(path, input, force)?;
    }
    Ok(())
}

fn check_target(path: &Path, input: &Path, force: bool) -> Result<Option<fs::Permissions>> {
    match fs::symlink_metadata(path) {
        Ok(metadata) => {
            if !force {
                bail!(
                    "refusing existing output target (no overwrite; use --force): {}",
                    path.display()
                );
            }
            ensure!(
                metadata.is_file(),
                "refusing non-regular output target (including directories and symlinks): {}",
                path.display()
            );
            ensure!(
                !same_file::is_same_file(input, path)
                    .with_context(|| format!("checking input alias {}", path.display()))?,
                "refusing to overwrite input through output target: {}",
                path.display()
            );
            Ok(Some(metadata.permissions()))
        }
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(e).with_context(|| format!("checking output {}", path.display())),
    }
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

    pub fn persist(self, input: &Path, force: bool) -> Result<PathBuf> {
        // Recheck target types and input aliases immediately before publication.
        // Rename replaces the directory entry, never writes through a target symlink.
        let finalize = || -> Result<()> {
            let permissions = check_target(&self.path, input, force)?;
            if force {
                if let Some(permissions) = permissions {
                    self.file
                        .as_file()
                        .set_permissions(permissions)
                        .with_context(|| {
                            format!("preserving output permissions for {}", self.path.display())
                        })?;
                }
                self.file.persist(&self.path)?;
            } else {
                self.file.persist_noclobber(&self.path)?;
            }
            Ok(())
        };
        finalize().with_context(|| {
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
        preflight(std::slice::from_ref(&target), &target, false).unwrap();
        let mut staged = Staged::new(target.clone()).unwrap();
        staged.file.write_all(b"our alignment").unwrap();
        fs::write(&target, b"racing creator").unwrap();
        assert!(staged.persist(&target, false).is_err());
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

    #[test]
    fn force_partial_publication_never_rolls_back_or_deletes_user_outputs() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("input");
        fs::write(&input, b"input").unwrap();
        let targets: Vec<_> = ["first", "second", "third"]
            .map(|name| dir.path().join(name))
            .into();
        for path in &targets {
            fs::write(path, b"original").unwrap();
        }
        preflight(&targets, &input, true).unwrap();
        let staged: Vec<_> = targets
            .iter()
            .map(|path| {
                let mut file = Staged::new(path.clone()).unwrap();
                file.file.write_all(b"replacement").unwrap();
                file
            })
            .collect();
        // Simulate a late filesystem change after staging. Publication of the
        // first file succeeds, the second fails, and the third is never attempted.
        fs::remove_file(&targets[1]).unwrap();
        fs::create_dir(&targets[1]).unwrap();
        let result: Result<Vec<_>> = staged
            .into_iter()
            .map(|file| file.persist(&input, true))
            .collect();
        let error = result.unwrap_err().to_string();
        assert!(error.contains("earlier outputs may already have been finalized"));
        assert_eq!(fs::read(&targets[0]).unwrap(), b"replacement");
        assert!(targets[1].is_dir());
        assert_eq!(fs::read(&targets[2]).unwrap(), b"original");
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 4);
    }

    #[test]
    fn failed_force_persist_retains_existing_target() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("input");
        let target = dir.path().join("alignment.phy");
        fs::write(&input, b"input").unwrap();
        fs::write(&target, b"original").unwrap();
        let staged = Staged::new(target.clone()).unwrap();
        // Force rename itself to fail, rather than the preflight check.
        fs::remove_file(staged.file.path()).unwrap();
        assert!(staged.persist(&input, true).is_err());
        assert_eq!(fs::read(target).unwrap(), b"original");
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 2);
    }
}
