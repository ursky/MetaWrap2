"""binning --semibin2: the command it builds, and the bins it hands over.

Checked on the command string, since both failures complete without error. ``multi_easy_bin``
takes the same inputs as ``single_easy_bin`` but expects a sample prefix on every contig name,
which MetaWrap2's assemblies lack. SemiBin2 also gzips its bins by default, so ``renumber_bins``
finds no FASTA and the bin set is empty.
"""

import os

from metawrap2.constants import bin_filename
from metawrap2.modules import binning


class _Args:
    def __init__(self, output, threads=4, min_len=1000, device=binning.DEFAULT_LEARNED_DEVICE):
        self.output = output
        self.threads = threads
        self.min_len = min_len
        self.device = device


def _capture(monkeypatch):
    """Replace the module's run() with a recorder, returning the list it appends to."""
    calls = []

    def fake_run(cmd, **kw):
        calls.append((cmd, kw))

    monkeypatch.setattr(binning, "run", fake_run)
    return calls


def _work_with_bins(tmp_path, names):
    work = tmp_path / "work_files"
    produced = work / "semibin2_out" / binning.SEMIBIN2_BIN_DIR
    produced.mkdir(parents=True)
    for name in names:
        (produced / name).write_text(">c1\nACGT\n")
    (work / "sample_0.bam").write_text("")
    return work


def test_the_command_is_single_easy_bin_and_does_not_gzip_its_bins(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    work = _work_with_bins(tmp_path, ["bin.1.fa", "bin.2.fa", "bin.10.fa"])
    out = tmp_path / "out"

    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">long\n" + "ACGT" * 1000 + "\n")

    binning._run_semibin2(
        _Args(str(out)),
        "e",
        str(work),
        str(assembly),
    )

    cmd = calls[0][0]
    assert "single_easy_bin" in cmd and "multi_easy_bin" not in cmd
    # Without this the bins are .fa.gz and the move below finds no FASTA.
    assert "--compression none" in cmd
    # Fixed seed: SemiBin2 trains a network, so otherwise two runs on one input give different bins
    # and benchmarks include that noise.
    assert "--random-seed %d" % binning.SEMIBIN2_SEED in cmd
    # CPU by default, for reproducibility and because defaults are costed in CPU-hours; otherwise
    # SemiBin2 trains on any GPU torch finds.
    assert "--engine cpu" in cmd and binning.DEFAULT_LEARNED_DEVICE == "cpu"

    # Bins are moved and renamed as for other binners (bin.10 after bin.2), leaving nothing in
    # SemiBin2's tree to be scored twice.
    assert sorted(os.listdir(out / "semibin2_bins")) == [bin_filename(i) for i in (1, 2, 3)]
    assert os.listdir(work / "semibin2_out" / binning.SEMIBIN2_BIN_DIR) == []


def test_semibin2_skips_short_assembly_without_running(tmp_path, monkeypatch):
    """SemiBin2 should be skipped if no contigs reach 4000 bp."""
    calls = _capture(monkeypatch)

    work = _work_with_bins(tmp_path, [])
    out = tmp_path / "out"

    assembly = tmp_path / "short.fa"
    assembly.write_text(">contig1\n" + "A" * 3999 + "\n>contig2\n" + "C" * 1000 + "\n")

    binning._run_semibin2(
        _Args(str(out)),
        "e",
        str(work),
        str(assembly),
    )

    # SemiBin2 must never be invoked.
    assert calls == []

    # Its output directory exists but contains no bins.
    bins_dir = out / "semibin2_bins"
    assert bins_dir.is_dir()
    assert list(bins_dir.iterdir()) == []


def test_semibin2_accepts_exactly_4000bp_contig(tmp_path, monkeypatch):
    """A contig exactly 4000 bp long passes the preflight."""
    calls = _capture(monkeypatch)

    work = _work_with_bins(tmp_path, ["bin.1.fa"])
    out = tmp_path / "out"

    assembly = tmp_path / "eligible.fa"
    assembly.write_text(">contig1\n" + "A" * 2000 + "\n" + "C" * 2000 + "\n")

    binning._run_semibin2(
        _Args(str(out)),
        "e",
        str(work),
        str(assembly),
    )

    assert len(calls) == 1
    assert "single_easy_bin" in calls[0][0]

    assert (out / "semibin2_bins" / bin_filename(1)).exists()


def test_semibin2_dry_run_does_not_require_assembly(tmp_path, monkeypatch):
    """Dry-run should work even when the staged FASTA does not exist."""
    calls = _capture(monkeypatch)

    work = _work_with_bins(tmp_path, [])

    monkeypatch.setattr(binning, "dry_run", lambda: True)

    binning._run_semibin2(
        _Args(str(tmp_path / "out")),
        "e",
        str(work),
        "not-yet-copied.fa",
    )

    assert len(calls) == 1
