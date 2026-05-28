"""Unit tests for individual transcript_tools functions.

These complement the end-to-end golden tests (test_transcript_tools_golden.py)
by exercising internal helpers directly, including edge cases that the example
fixtures do not trigger.
"""

import types

from leafcutter2 import transcript_tools as tt


# ── A5: Is_bedline_complete robustness ─────────────────────────────────────

def test_is_bedline_complete_handles_trailing_comma():
    """Some bedparse versions emit a trailing comma; int('') must not be hit."""
    bl = types.SimpleNamespace(exStarts="0,100,", exLengths="50,50,", start=1000, end=1150)
    assert tt.Is_bedline_complete(bl) is True


def test_is_bedline_complete_without_trailing_comma():
    bl = types.SimpleNamespace(exStarts="0,100", exLengths="50,50", start=1000, end=1150)
    assert tt.Is_bedline_complete(bl) is True


def test_is_bedline_complete_detects_incomplete():
    # last block ends before transcript end -> incomplete
    bl = types.SimpleNamespace(exStarts="0,100", exLengths="50,50", start=1000, end=2000)
    assert tt.Is_bedline_complete(bl) is False


# ── A4 / NMDFinderB classification branches ────────────────────────────────
# Sequences are marked with '^' start, '*' stop, '|' splice junction.
# Expected values verified against the classifier's branch logic.

def test_nmd_no_cds():
    assert tt.get_NMD_detective_B_classification("ACGTACG|CACGT") == "No CDS"


def test_nmd_no_stop():
    assert tt.get_NMD_detective_B_classification("A^ATGACG|CACGT") == "No stop"


def test_nmd_last_exon():
    # stop in the final exon
    assert tt.get_NMD_detective_B_classification("ACGT|A^CGT*ACG") == "Last exon"


def test_nmd_start_proximal():
    # short CDS (<=125 nt) with an internal (non-final) stop exon
    assert tt.get_NMD_detective_B_classification("ACG^TAAAAAA*AAA|CACGT") == "Start proximal"


def test_nmd_long_exon():
    # CDS > 125 nt and the stop-containing exon >= 407 nt
    seq = "ACG^T" + "A" * 200 + "|" + "A" * 407 + "*A|CACGT"
    assert tt.get_NMD_detective_B_classification(seq) == "Long exon"


def test_nmd_50nt_rule():
    # stop within 50 nt of the last junction
    seq = "ACG^T" + "A" * 200 + "*" + "A" * 30 + "|CACGT"
    assert tt.get_NMD_detective_B_classification(seq) == "50 nt rule"


def test_nmd_trigger():
    # internal stop, far from last junction, exon < 407
    seq = "ACG^T" + "A" * 200 + "*" + "A" * 200 + "|CACGT"
    assert tt.get_NMD_detective_B_classification(seq) == "Trigger NMD"


def test_nmd_does_not_crash_on_guarded_branch():
    # A4: previously the InternalStopExon match could be None at the "Long exon"
    # branch and crash. Must return a string, not raise.
    seq = "ACG^T" + "A" * 200 + "*" + "A" * 60 + "|CACGT"
    assert isinstance(tt.get_NMD_detective_B_classification(seq), str)
