"""Run with python3 -m unittest discover -s scripts -p 'test_*.py'."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import benchmark


class BenchmarkTests(unittest.TestCase):
    def test_generator_deterministic(self):
        with tempfile.TemporaryDirectory() as folder:
            first, second = Path(folder) / "a", Path(folder) / "b"
            expected = benchmark.generate(first, 3, 7, 42)
            self.assertEqual(expected, benchmark.generate(second, 3, 7, 42))
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(len(expected[0]), 3)

    def test_disk_usage(self):
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder) / "test").write_bytes(b"123")
            self.assertEqual(benchmark.disk_usage(folder)[0], 3)

    def test_sampling_failure_is_raised_after_child_reaped(self):
        with tempfile.TemporaryDirectory() as folder:
            run_dir = Path(folder)
            (run_dir / "outputs").mkdir()
            (run_dir / "spool").mkdir()
            error = OSError("injected sampling failure")
            processes = []

            def start_process(*args, **kwargs):
                process = real_popen(*args, **kwargs)
                processes.append(process)
                return process

            real_popen = subprocess.Popen
            with mock.patch.object(benchmark, "disk_usage", side_effect=error), \
                 mock.patch.object(benchmark.subprocess, "Popen", side_effect=start_process):
                with self.assertRaisesRegex(RuntimeError, "disk sampling failed") as raised:
                    benchmark.measure([sys.executable, "-c", "pass"], run_dir, 0.01)
            self.assertIs(raised.exception.__cause__, error)
            self.assertEqual(processes[0].returncode, 0)

    def test_rust_smoke_and_detect_corruption(self):
        binary = Path(os.environ.get("VCF2PHYLIP_BENCH_BINARY", benchmark.ROOT / "target/release/vcf2phylip-rs"))
        if not binary.exists():
            self.skipTest("build release binary first or set VCF2PHYLIP_BENCH_BINARY")
        with tempfile.TemporaryDirectory(dir=benchmark.ROOT / "target") as folder:
            report = benchmark.main(["--rust-tool", str(binary), "--sandbox-parent", folder,
                                     "--samples", "3", "--sites", "83", "--memory-limit", "6"])
            result = json.loads(report.read_text())
            self.assertTrue(result["validated"])
            self.assertIsNone(result["baseline"])
            self.assertEqual([r["mode"] for r in result["runs"]], ["memory", "disk"])
            expected, sites = benchmark.generate(Path(folder) / "expected", 3, 83, 42)
            output = report.parent / "disk-0/outputs"
            (output / "synthetic.min1.fasta").write_text(">S0\nBAD\n")
            with self.assertRaisesRegex(ValueError, "FASTA mismatch"):
                benchmark.validate_outputs(output, expected, sites)


if __name__ == "__main__":
    unittest.main()
