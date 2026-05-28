"""Golden-file regression tests for leafcutter2-transcript-tools.

Each matrix row runs the CLI on a small committed GTF fixture across a
combination of settings, then compares the produced GTF and BED12 against
committed golden files in tests/golden/transcript_tools/.

Comparison is on decompressed, comment-stripped content: the GTF output
carries a "#! args: ..." header containing absolute paths that vary per run
and per machine, so all '#'-prefixed lines are filtered before comparing.

To regenerate goldens after a reviewed/approved change:

    mamba run -n leafcutter2 python tests/test_transcript_tools_golden.py

Only regenerate when you have inspected the diff and confirmed the change is
intentional.
"""

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent
FIXTURE_DIR = Path(__file__).parent / "fixtures" / "transcript_tools"
GOLDEN_DIR = Path(__file__).parent / "golden" / "transcript_tools"
FASTA = REPO_ROOT / "example" / "annotation" / "chr10.fa.gz"

WCDS = "chr10_minimal_w_CDS.subset.gtf"
NOCDS = "chr10_minimal.subset.gtf"

# (key, fixture, extra CLI args). Each row emits one GTF + one BED12 golden.
MATRIX = [
    ("wCDS_approachA", WCDS, ["-translation_approach", "A"]),
    ("wCDS_approachB", WCDS, ["-translation_approach", "B"]),
    ("wCDS_approachC", WCDS, ["-translation_approach", "C"]),
    ("wCDS_approachE", WCDS, ["-translation_approach", "E"]),
    ("wCDS_approachF", WCDS, ["-translation_approach", "F"]),
    ("noCDS_approachD", NOCDS, ["-translation_approach", "D"]),
    ("wCDS_approachB_uorf", WCDS, ["-translation_approach", "B", "--include_uorf_analysis"]),
    ("wCDS_approachB_sorted", WCDS, ["-translation_approach", "B", "--sort_output"]),
]


def _strip_comments(text: str) -> str:
    return "\n".join(line for line in text.splitlines() if not line.startswith("#"))


def run_matrix_row(key: str, fixture: str, extra_args: list, out_dir: Path):
    """Run the CLI for one matrix row; return (gtf_content, bed12_content), comment-stripped."""
    gtf_out = out_dir / f"{key}.gtf"
    bed12_out = out_dir / f"{key}.bed12"
    cmd = [
        "leafcutter2-transcript-tools",
        "-i", str(FIXTURE_DIR / fixture),
        "-fa", str(FASTA),
        "-o", str(gtf_out),
        "-bed12_out", str(bed12_out),
        *extra_args,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, (
        f"transcript-tools failed for {key}.\ncmd: {' '.join(cmd)}\n"
        f"stdout: {result.stdout}\nstderr: {result.stderr}"
    )
    return _strip_comments(gtf_out.read_text()), _strip_comments(bed12_out.read_text())


@pytest.mark.skipif(not FASTA.exists(), reason="Example genome FASTA not found")
@pytest.mark.parametrize("key,fixture,extra_args", MATRIX, ids=[m[0] for m in MATRIX])
def test_transcript_tools_golden(key, fixture, extra_args, tmp_path):
    gtf_golden = GOLDEN_DIR / f"{key}.gtf"
    bed12_golden = GOLDEN_DIR / f"{key}.bed12"
    if not gtf_golden.exists() or not bed12_golden.exists():
        pytest.skip(f"golden missing for {key}; run this file as a script to generate")

    gtf_content, bed12_content = run_matrix_row(key, fixture, extra_args, tmp_path)

    assert gtf_content == gtf_golden.read_text(), f"GTF output differs from golden for {key}"
    assert bed12_content == bed12_golden.read_text(), f"BED12 output differs from golden for {key}"


def _regenerate():
    import tempfile

    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        out_dir = Path(td)
        for key, fixture, extra_args in MATRIX:
            gtf_content, bed12_content = run_matrix_row(key, fixture, extra_args, out_dir)
            (GOLDEN_DIR / f"{key}.gtf").write_text(gtf_content)
            (GOLDEN_DIR / f"{key}.bed12").write_text(bed12_content)
            print(f"wrote golden: {key}")


if __name__ == "__main__":
    _regenerate()
