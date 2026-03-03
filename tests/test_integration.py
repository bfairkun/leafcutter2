"""Integration tests for leafcutter2.

These tests run the full pipeline (or sub-steps) against the example data
bundled in the example/ directory of the repository.
"""

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

# Paths to example data (relative to repo root)
REPO_ROOT = Path(__file__).parent.parent
EXAMPLE_DIR = REPO_ROOT / "example"
ANNOTATION_GTF = EXAMPLE_DIR / "annotation" / "chr10.gtf.gz"
ANNOTATION_FA = EXAMPLE_DIR / "annotation" / "chr10.fa.gz"
JUNCTION_FILES_TXT = EXAMPLE_DIR / "junction_files.txt"


def _make_abs_juncfiles(tmp_path: Path) -> Path:
    """Write a junction_files.txt with absolute paths into tmp_path."""
    abs_junc_txt = tmp_path / "junction_files.txt"
    lines = []
    for line in JUNCTION_FILES_TXT.read_text().splitlines():
        line = line.strip()
        if line:
            lines.append(str(EXAMPLE_DIR / line))
    abs_junc_txt.write_text("\n".join(lines) + "\n")
    return abs_junc_txt


def test_package_imports():
    """The package and all sub-modules import without error."""
    import importlib
    for module in [
        "leafcutter2",
        "leafcutter2.cli",
        "leafcutter2.star_utils",
    ]:
        importlib.import_module(module)


def test_entry_points_help():
    """All registered entry-point commands respond to --help."""
    for cmd in ["leafcutter2", "leafcutter2-make-clusters", "leafcutter2-star2junc"]:
        result = subprocess.run(
            [cmd, "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"`{cmd} --help` exited with code {result.returncode}.\n"
            f"stderr: {result.stderr}"
        )


def test_make_clusters(tmp_path):
    """leafcutter2-make-clusters runs on the example data and produces a clusters file."""
    abs_junc_txt = _make_abs_juncfiles(tmp_path)
    result = subprocess.run(
        [
            "leafcutter2-make-clusters",
            "-j", str(abs_junc_txt),
            "-r", str(tmp_path),
            "-o", "leafcutter2",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"leafcutter2-make-clusters failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    clusters_file = tmp_path / "clustering" / "leafcutter2_clusters"
    assert clusters_file.exists(), f"Expected clusters file not found: {clusters_file}"
    assert clusters_file.stat().st_size > 0, "Clusters file is empty"


@pytest.mark.skipif(
    not ANNOTATION_GTF.exists() or not ANNOTATION_FA.exists(),
    reason="Example annotation files not found",
)
def test_full_pipeline(tmp_path):
    """Full leafcutter2 pipeline runs on example data and produces output files."""
    abs_junc_txt = _make_abs_juncfiles(tmp_path)
    result = subprocess.run(
        [
            "leafcutter2",
            "-j", str(abs_junc_txt),
            "-r", str(tmp_path),
            "-o", "leafcutter2",
            "-A", str(ANNOTATION_GTF),
            "-G", str(ANNOTATION_FA),
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"leafcutter2 pipeline failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    expected_files = [
        tmp_path / "leafcutter2.cluster_ratios.gz",
        tmp_path / "leafcutter2.junction_counts.gz",
    ]
    for f in expected_files:
        assert f.exists(), f"Expected output file not found: {f}"
        assert f.stat().st_size > 0, f"Output file is empty: {f}"

    # Check that output files are valid gzip with at least a header line
    for f in expected_files:
        with gzip.open(f, "rt") as fh:
            header = fh.readline()
        assert header.strip(), f"Output file has no header: {f}"
