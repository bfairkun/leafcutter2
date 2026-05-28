"""Unit tests for individual transcript_tools functions.

These complement the end-to-end golden tests (test_transcript_tools_golden.py)
by exercising internal helpers directly, including edge cases that the example
fixtures do not trigger.
"""

import types

import bedparse

from leafcutter2 import transcript_tools as tt


def _make_tx(strand, cds_start=0, cds_end=250):
    """3 exons x 50bp with 50bp introns; CDS span configurable.

    BED fields: chr start end name score strand cdsStart cdsEnd color nEx
    exLengths exStarts. Exons span genomic 0-50, 100-150, 200-250 (tx length 150).
    """
    return bedparse.bedline(
        ["chr1", 0, 250, "tx", 0, strand, cds_start, cds_end, 0, 3, "50,50,50", "0,100,200"]
    )


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


# ── C5: calculate_frames (GTF phase) for both strands ──────────────────────

def test_calculate_frames_plus_strand():
    # 3 CDS blocks of 50nt each; phase of the first translated block is 0.
    assert tt.calculate_frames(_make_tx("+")) == [0, 1, 2]


def test_calculate_frames_minus_strand():
    # Genomic order; the 3'-most genomic block is the translation start (phase 0).
    assert tt.calculate_frames(_make_tx("-")) == [2, 1, 0]


def test_calculate_frames_translation_start_phase_zero():
    plus = tt.calculate_frames(_make_tx("+"))
    minus = tt.calculate_frames(_make_tx("-"))
    assert plus[0] == 0          # + strand: first genomic block starts translation
    assert minus[-1] == 0        # - strand: last genomic block starts translation
    assert all(p in (0, 1, 2) for p in plus + minus)


# ── C5: extract_codon places 3bp at the correct genomic end per strand ─────

def test_extract_codon_plus_strand():
    bl = _make_tx("+")
    start = tt.extract_codon(bl, "start")
    stop = tt.extract_codon(bl, "stop")
    assert (start.start, start.end) == (0, 3)        # 5' genomic end
    assert (stop.start, stop.end) == (247, 250)      # 3' genomic end


def test_extract_codon_minus_strand():
    bl = _make_tx("-")
    start = tt.extract_codon(bl, "start")
    stop = tt.extract_codon(bl, "stop")
    assert (start.start, start.end) == (247, 250)    # start is at higher genomic coord
    assert (stop.start, stop.end) == (0, 3)


def test_extract_codon_requires_cds():
    # transcript with no CDS (cdsStart == cdsEnd) -> ValueError
    import pytest
    bl = bedparse.bedline(["chr1", 0, 250, "tx", 0, "+", 0, 0, 0, 3, "50,50,50", "0,100,200"])
    with pytest.raises(ValueError):
        tt.extract_codon(bl, "start")


# ── C5: get_absolute_pos transcript->genomic mapping ───────────────────────

def test_get_absolute_pos_plus_strand():
    bl = _make_tx("+")
    assert tt.get_absolute_pos(bl, 0) == 0          # first tx base -> 5' genomic
    assert tt.get_absolute_pos(bl, 149) == 249      # last tx base -> 3' genomic


def test_get_absolute_pos_minus_strand_monotonic_decreasing():
    # On the minus strand, increasing transcript coordinate maps to decreasing
    # genomic coordinate.
    bl = _make_tx("-")
    assert tt.get_absolute_pos(bl, 0) > tt.get_absolute_pos(bl, 149)
