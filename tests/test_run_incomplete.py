"""RUN_INCOMPLETE: a run that stops early says so at the top of its output folder.

A run stopped by an error, the disk ceiling or a signal leaves a partial tree, and a tool that
globs ``metawrap_*_bins`` cannot tell it from a run that recovered nothing. The marker is written
when the run starts, rewritten with the step and the reason if it stops early, and removed when it
finishes. bin_refinement is the module under test; the marker is written by the shared
``start_run``/``finish_run``, so every module gets it.
"""

from __future__ import annotations

import os

import pytest
import test_consolidator as helpers

from metawrap2 import checkm, cli, limits
from metawrap2.modules import _common

MARKER = _common.INCOMPLETE_MARKER


@pytest.fixture
def fake_checkm(scored_dirs):
    return scored_dirs


scored_dirs = helpers.scored_dirs


def _argv(root, extra=()):
    dirs = helpers._inputs(root)
    out = root / "refined"
    argv = ["bin_refinement", "-o", str(out), "-t", "1", "-m", "8", "-x", "10", "--skip-gunc"]
    return argv + ["-A", dirs[0], "-B", dirs[1], *extra], str(out)


def _marker(out):
    with open(os.path.join(out, MARKER)) as fh:
        return fh.read()


def test_a_finished_run_leaves_no_marker(tmp_path, fake_checkm):
    argv, out = _argv(tmp_path, ["-c", "50"])
    assert cli.main(argv) == 0
    assert not os.path.exists(os.path.join(out, MARKER))
    assert os.path.isdir(os.path.join(out, "metawrap_50_10_bins"))


def test_an_error_leaves_the_step_and_the_reason_and_exits_non_zero(tmp_path, fake_checkm, capsys):
    # No bin is more than 100% complete, so the run stops in its finalize step.
    argv, out = _argv(tmp_path, ["-c", "100"])
    with pytest.raises(SystemExit) as stopped:
        cli.main(argv)
    assert stopped.value.code == 1
    text = _marker(out)
    assert "metawrap2 bin_refinement stopped before finishing (failed)" in text
    assert "Step: finalize" in text and "Last finished step: consolidate" in text
    assert "Reason: Consolidation produced no bins" in text
    assert "--resume" in text
    err = capsys.readouterr().err
    assert "metawrap2 bin_refinement: error: Consolidation produced no bins" in err


def test_a_disk_stop_is_reported_and_a_resume_removes_the_marker(
    tmp_path, monkeypatch, fake_checkm, capsys
):
    real = checkm.run_checkm

    def disk_fills(bins_dir, **kwargs):
        monitor = limits.DiskMonitor([], limits.resolve())
        monitor.tripped = limits.DiskState(
            path=str(tmp_path), total_bytes=1000, used_bytes=950, ceiling_bytes=700
        )
        monkeypatch.setattr(_common, "_active_disk_monitor", monitor)
        raise KeyboardInterrupt

    monkeypatch.setattr(checkm, "run_checkm", disk_fills)
    argv, out = _argv(tmp_path, ["-c", "50"])
    assert cli.main(argv) == 130
    text = _marker(out)
    assert "stopped before finishing (interrupted)" in text
    assert "Step: checkm" in text
    assert "Reason: the disk ceiling was reached:" in text and "95% full" in text
    assert "Stopped before finishing" in capsys.readouterr().err
    assert not os.path.exists(os.path.join(out, "metawrap_50_10_bins"))

    monkeypatch.setattr(checkm, "run_checkm", real)
    assert cli.main(["--resume"] + argv) == 0
    assert not os.path.exists(os.path.join(out, MARKER))
    assert os.path.isdir(os.path.join(out, "metawrap_50_10_bins"))


def test_the_marker_alone_does_not_count_as_an_earlier_run(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / MARKER).write_text("metawrap2 x started\n")
    assert not _common._has_prior_output(str(out))
