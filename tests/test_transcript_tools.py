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
