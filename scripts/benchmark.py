#!/usr/bin/env python3
"""Opt-in, stdlib-only synthetic benchmark. Never downloads a baseline or builds."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import shutil
import subprocess
import sys
import tempfile
import threading
import time

BASELINE_COMMIT = "85a260e655b5f3c5dc2cb9915eb24552a002357d"
ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_text(argv):
    try:
        return subprocess.check_output(argv, cwd=ROOT, stderr=subprocess.STDOUT,
                                       text=True).strip()
    except (OSError, subprocess.CalledProcessError) as error:
        return f"unavailable: {error}"


def hardware():
    result = dict(machine=platform.machine(), processor=platform.processor(),
                  cpu_count=os.cpu_count())
    if sys.platform == "linux":
        for source, key, field in [("/proc/cpuinfo", "cpu_model", "model name"),
                                   ("/proc/meminfo", "host_mem_total", "MemTotal")]:
            try:
                with open(source) as stream:
                    for line in stream:
                        name, separator, value = line.partition(":")
                        if separator and name.strip() == field:
                            result[key] = value.strip()
                            break
            except OSError:
                result[key] = "unavailable"
    return result


def generate(path, samples, sites, seed):
    """Canonical biallelic SNPs, fully called diploids or ./.; sample0 called."""
    rng = random.Random(seed)
    digests = [hashlib.sha256() for _ in range(samples)]
    site_digest = hashlib.sha256()
    with path.open("w") as stream:
        stream.write('##fileformat=VCFv4.2\n##contig=<ID=chr1>\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                     "\t".join(f"S{i}" for i in range(samples)) + "\n")
        for site in range(1, sites + 1):
            calls = ["0/0"] + [rng.choice(["0/0", "0/1", "1/1", "./."])
                                 for _ in range(samples - 1)]
            symbols = {"0/0": b"A", "0/1": b"M", "1/1": b"C", "./.": b"N"}
            for digest, call in zip(digests, calls):
                digest.update(symbols[call])
            called = sum(call != "./." for call in calls)
            site_digest.update(f"chr1\t{site}\t{called}\n".encode())
            stream.write(f"chr1\t{site}\t.\tA\tC\t.\tPASS\t.\tGT\t" +
                         "\t".join(calls) + "\n")
    return [(f"S{i}", sites, digest.hexdigest()) for i, digest in enumerate(digests)], site_digest.hexdigest()


def disk_usage(path):
    logical = allocated = 0
    def raise_walk_error(error):
        raise error

    for root, _, files in os.walk(path, onerror=raise_walk_error):
        for name in files:
            try:
                stat = (Path(root) / name).stat()
            except FileNotFoundError:
                continue  # rename/unlink between directory listing and stat
            logical += stat.st_size
            allocated += getattr(stat, "st_blocks", 0) * 512
    return logical, allocated


def measure(argv, run_dir, interval):
    output, spool = run_dir / "outputs", run_dir / "spool"
    peaks = {"output_and_staging": [0, 0], "matrix_spool": [0, 0]}
    stop = threading.Event()
    sampling_errors = []

    def sample():
        for label, path in [("output_and_staging", output), ("matrix_spool", spool)]:
            usage = disk_usage(path)
            peaks[label] = [max(a, b) for a, b in zip(peaks[label], usage)]

    def monitor():
        try:
            while not stop.is_set():
                sample()
                stop.wait(interval)
            sample()
        except Exception as error:
            sampling_errors.append(error)

    env = dict(os.environ, TMPDIR=str(spool), TMP=str(spool), TEMP=str(spool))
    with (run_dir / "stdout.txt").open("w") as stdout, (run_dir / "stderr.txt").open("w") as stderr:
        start = time.perf_counter()
        process = subprocess.Popen(argv, stdout=stdout, stderr=stderr, env=env, cwd=run_dir)
        thread = threading.Thread(target=monitor)
        thread.start()
        try:
            if hasattr(os, "wait4") and sys.platform in ("linux", "darwin"):
                _, status, usage = os.wait4(process.pid, 0)
                process.returncode = os.waitstatus_to_exitcode(status)
                rss = usage.ru_maxrss * (1024 if sys.platform == "linux" else 1)
                rss_method = "wait4 per-child ru_maxrss; Linux KiB/macOS bytes normalized to bytes"
            else:
                process.wait()
                rss, rss_method = None, "unsupported platform: no per-child RSS measurement"
            wall = time.perf_counter() - start
        finally:
            stop.set()
            thread.join()
    if sampling_errors:
        raise RuntimeError(f"disk sampling failed for command {argv} in {run_dir}") from sampling_errors[0]
    if process.returncode:
        raise RuntimeError(f"command failed ({process.returncode}): {argv}; see {run_dir}/stderr.txt")
    return dict(argv=argv, wall_seconds=wall, peak_rss_bytes=rss, rss_method=rss_method,
                sampled_disk_peaks_bytes=peaks, final_output_bytes=disk_usage(output),
                stdout=(run_dir / "stdout.txt").read_text())


def validate_outputs(folder, expected, expected_sites):
    paths = list(folder.iterdir())
    fasta = next(p for p in paths if p.suffix == ".fasta")
    phy = next(p for p in paths if p.suffix == ".phy")
    used = next(p for p in paths if p.name.endswith("used_sites.tsv"))
    rows, name, digest, length = [], None, None, 0
    with fasta.open() as stream:
        for line in stream:
            if line.startswith(">"):
                if name is not None:
                    rows.append((name, length, digest.hexdigest()))
                name, digest, length = line[1:].strip(), hashlib.sha256(), 0
            else:
                sequence = "".join(line.split()).upper().encode()
                if digest is None:
                    raise ValueError("FASTA sequence before identifier")
                digest.update(sequence)
                length += len(sequence)
    if name is not None:
        rows.append((name, length, digest.hexdigest()))
    if rows != expected:
        raise ValueError(f"FASTA mismatch: {fasta}")
    with phy.open() as stream:
        dimensions = tuple(map(int, next(stream).split()))
        if dimensions != (len(expected), expected[0][1]):
            raise ValueError("PHYLIP dimensions mismatch")
        rows = []
        for line in stream:
            label, sequence = line.split(maxsplit=1)
            sequence = "".join(sequence.split()).upper().encode()
            rows.append((label, len(sequence), hashlib.sha256(sequence).hexdigest()))
        if rows != expected:
            raise ValueError(f"PHYLIP mismatch: {phy}")
    digest = hashlib.sha256()
    with used.open() as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError("invalid site report")
            digest.update(("\t".join(fields) + "\n").encode())
    if digest.hexdigest() != expected_sites:
        raise ValueError(f"retained-site mismatch: {used}")
    return {p.name: sha256(p) for p in paths if p.is_file()}


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rust-tool", type=Path, default=ROOT / "target/release/vcf2phylip-rs")
    parser.add_argument("--python-tool", type=Path, help="local upstream v2.9 script; never fetched")
    parser.add_argument("--python-executable", default=sys.executable)
    parser.add_argument("--samples", type=positive, default=8)
    parser.add_argument("--sites", type=positive, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--repeats", type=positive, default=1)
    parser.add_argument("--memory-limit", default="64KiB")
    parser.add_argument("--poll-seconds", type=float, default=0.01)
    parser.add_argument("--sandbox-parent", type=Path, default=ROOT / "target/bench")
    args = parser.parse_args(argv)
    if args.poll_seconds <= 0 or not args.poll_seconds < float("inf"):
        parser.error("--poll-seconds must be finite and positive")
    rust = args.rust_tool.resolve(strict=True)
    python_tool = args.python_tool.resolve(strict=True) if args.python_tool else None
    python_executable = shutil.which(args.python_executable)
    if python_tool and not python_executable:
        parser.error("Python executable not found")
    args.sandbox_parent.mkdir(parents=True, exist_ok=True)
    sandbox = Path(tempfile.mkdtemp(prefix="bench-", dir=args.sandbox_parent.resolve()))
    source = sandbox / "synthetic.vcf"
    expected, expected_sites = generate(source, args.samples, args.sites, args.seed)
    metadata = dict(settings=vars(args), sandbox=sandbox, input_sha256=sha256(source),
                    source_commit=command_text(["git", "rev-parse", "HEAD"]),
                    source_dirty=command_text(["git", "status", "--porcelain"]),
                    lock_sha256=sha256(ROOT / "Cargo.lock"), rust_binary_sha256=sha256(rust),
                    rust_version=command_text(["rustc", "--version"]),
                    cargo_version=command_text(["cargo", "--version"]),
                    tool_version=command_text([str(rust), "--version"]),
                    python_version=sys.version, platform=platform.platform(),
                    hardware=hardware(), uname=list(platform.uname()),
                    baseline=None if python_tool is None else dict(path=python_tool,
                        sha256=sha256(python_tool), expected_commit=BASELINE_COMMIT,
                        provenance="user-supplied; commit identity not verified by file hash alone",
                        interpreter=command_text([python_executable, "--version"])),
                    comparison="rust-only smoke; no Python speedup claim" if python_tool is None else
                               "local supplied baseline; comparable only after all validations pass",
                    disk_method="polling independent logical/allocated peaks; may miss short-lived files; output peak includes final files; final outputs also reported separately; peaks not simultaneous",
                    output_dimensions=[args.samples, args.sites],
                    rust_matrix_payload_bytes=args.samples * args.sites,
                    validated=False, runs=[])
    report_path = sandbox / "report.json"
    try:
        for repeat in range(args.repeats):
            for mode in ["memory", "disk"] + (["python"] if python_tool else []):
                run_dir = sandbox / f"{mode}-{repeat}"
                output, spool = run_dir / "outputs", run_dir / "spool"
                output.mkdir(parents=True)
                spool.mkdir()
                common = ["-i", str(source), "-m", "1", "-f", "-w", "--output-folder", str(output), "--output-prefix", "synthetic"]
                if mode == "python":
                    command = [python_executable, str(python_tool)] + common
                else:
                    # Match the baseline's count-only filtering and .min1 names.
                    command = [str(rust)] + common + ["--min-call-rate", "0", "--storage", mode]
                    if mode == "disk":
                        command += ["--memory-limit", args.memory_limit, "--temp-dir", str(spool)]
                metadata["pending_run"] = dict(argv=command, mode=mode, repeat=repeat)
                result = measure(command, run_dir, args.poll_seconds)
                result.update(mode=mode, repeat=repeat)
                result["output_sha256"] = validate_outputs(output, expected, expected_sites)
                result["correctness"] = "FASTA/PHYLIP rows and dimensions + retained coordinates/counts match generator"
                metadata["runs"].append(result)
                del metadata["pending_run"]
        metadata["validated"] = True
    finally:
        report_path.write_text(json.dumps(metadata, indent=2, default=str) + "\n")
    with (sandbox / "runs.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["mode", "repeat", "wall_seconds", "peak_rss_bytes", "output_peak_logical_bytes", "output_peak_allocated_bytes", "spool_peak_logical_bytes", "spool_peak_allocated_bytes"])
        for run in metadata["runs"]:
            writer.writerow([run["mode"], run["repeat"], run["wall_seconds"], run["peak_rss_bytes"], *run["sampled_disk_peaks_bytes"]["output_and_staging"], *run["sampled_disk_peaks_bytes"]["matrix_spool"]])
    print(report_path)
    return report_path


if __name__ == "__main__":
    main()
