"""bin_refinement -b/--bin-sets: any number of bin sets, given comma-separated or repeated.

-A..-F still work (hidden from -h) but cannot be mixed with -b. The default consolidator and
dastool take any number of sets; --mw1-refinement takes at most MW1_MAX_BIN_SETS. CheckM is the
stand-in from test_consolidator, which scores a bin against the genome it mostly comes from.
"""

from __future__ import annotations

import os
import runpy
import shlex
import sys
import time

import pytest
import test_consolidator as helpers

from metawrap2 import command
from metawrap2.io.seqio import contig_id, iter_fasta
from metawrap2.modules import bin_refinement

scored_dirs = helpers.scored_dirs

#: Five binners over the two genomes of test_consolidator, each splitting them differently.
FIVE = [
    {"bin.1": ["c01", "c02", "c03", "c04"], "bin.2": ["c07", "c08", "c09", "c10"]},
    {"bin.1": ["c03", "c04", "c05", "c06"], "bin.2": ["c07", "c08", "c09", "c11"]},
    {"bin.1": ["c01", "c02", "c03"], "bin.2": ["c08", "c09", "c10"]},
    {"bin.1": ["c02", "c03", "c04", "c05"], "bin.2": ["c07", "c08", "c10", "c12"]},
    {"bin.1": ["c04", "c05", "c06"], "bin.2": ["c07", "c09", "c10"]},
]


def _sets(root, specs):
    dirs = []
    for index, spec in enumerate(specs):
        directory = root / ("binner_%d" % index)
        directory.mkdir(parents=True)
        for name, contigs in spec.items():
            with open(directory / ("%s.fa" % name), "w") as fh:
                fh.writelines(">%s\n%s\n" % (c, helpers.CONTIGS[c]) for c in contigs)
        dirs.append(str(directory))
    return dirs


def _collected(argv):
    args = bin_refinement._parse_args(["-o", "o"] + argv)
    return bin_refinement.collect_bin_sets(args)


def _refuses(argv, capsys):
    with pytest.raises(SystemExit):
        _collected(argv)
    captured = capsys.readouterr()
    return " ".join((captured.out + captured.err).replace("*", " ").split())


def test_the_comma_form_the_repeated_form_and_the_old_names_agree(tmp_path):
    dirs = _sets(tmp_path, FIVE[:3])
    expected = [os.path.abspath(d) for d in dirs]
    assert _collected(["-b", ",".join(dirs)]) == expected
    assert _collected(["-b", dirs[0], "-b", dirs[1], "--bin-sets", dirs[2]]) == expected
    assert _collected(["-b", ",".join(dirs[:2]), "-b", dirs[2]]) == expected
    assert _collected(["--bin_sets", ",".join(dirs)]) == expected
    assert _collected(["--bins", dirs[0], "--bins", dirs[1], "--bins", dirs[2]]) == expected
    assert _collected(["-A", dirs[0], "-B", dirs[1], "-C", dirs[2]]) == expected


def test_a_folder_with_a_comma_in_its_name_can_be_given_on_its_own(tmp_path):
    odd = tmp_path / "bins,metabat2"
    odd.mkdir()
    other = tmp_path / "semibin2"
    other.mkdir()
    assert _collected(["-b", str(odd), "-b", str(other)]) == [str(odd), str(other)]


def test_mixing_b_with_the_lettered_flags_is_refused(tmp_path, capsys):
    dirs = _sets(tmp_path, FIVE[:2])
    said = _refuses(["-A", dirs[0], "-b", dirs[1]], capsys)
    assert "Give the bin sets one way" in said and "-b %s,%s" % (dirs[0], dirs[1]) in said


def test_a_duplicate_or_a_missing_folder_is_refused(tmp_path, capsys, monkeypatch):
    dirs = _sets(tmp_path, FIVE[:2])
    said = _refuses(["-b", "%s,%s,%s/" % (dirs[0], dirs[1], dirs[0])], capsys)
    assert "is given twice" in said
    missing = str(tmp_path / "nope")
    said = _refuses(["-b", "%s,%s" % (dirs[0], missing)], capsys)
    assert "Not a folder: %s" % missing in said and "comma" in said
    assert "No bin sets given" in _refuses([], capsys)
    # A dry run plans folders that a study has not made yet.
    monkeypatch.setattr(command.runner, "dry_run", True)
    assert _collected(["-b", "%s,%s" % (dirs[0], missing)]) == [dirs[0], missing]


def test_the_lettered_flags_are_hidden_and_b_is_listed(capsys):
    with pytest.raises(SystemExit):
        bin_refinement._parse_args(["-h"])
    text = capsys.readouterr().out
    # Python 3.14 lists "-b, --bin-sets DIR1,DIR2,..."; earlier versions repeat the metavar.
    assert "--bin-sets DIR1,DIR2,..." in text and "-b" in text.split("--bin-sets", 1)[0]
    assert "--bins-a" not in text and "--bin_sets" not in text and "\n  --bins " not in text


def test_labels_continue_past_z():
    labels = [bin_refinement.bin_set_label(i) for i in (0, 1, 2, 25, 26, 27, 51, 52)]
    assert labels == [
        "binsA",
        "binsB",
        "binsC",
        "binsZ",
        "binsAA",
        "binsAB",
        "binsAZ",
        "binsBA",
    ]


@pytest.mark.parametrize("n_sets", [4, 5])
def test_the_default_consolidates_four_and_five_sets(tmp_path, scored_dirs, n_sets):
    dirs = _sets(tmp_path, FIVE[:n_sets])
    out = tmp_path / "refined"
    argv = ["-o", str(out), "-t", "1", "-m", "8", "-c", "50", "-x", "10", "--skip-gunc"]
    assert bin_refinement.main(argv + ["-b", ",".join(dirs)]) == 0
    final = helpers._final_bins(str(out))
    assert sorted(final.values()) == [helpers.GENOMES["g1"], helpers.GENOMES["g2"]]
    work = os.listdir(os.path.join(str(out), "work_files"))
    labels = [bin_refinement.bin_set_label(i) for i in range(n_sets)]
    assert set(labels) <= set(work)
    assert set(labels) <= set(scored_dirs)


def test_mw1_refinement_refuses_more_sets_than_it_takes(tmp_path, scored_dirs, capsys):
    specs = [{"bin.1": ["c%02d" % (i + 1), "c%02d" % (i + 2)]} for i in range(7)]
    dirs = _sets(tmp_path, specs)
    argv = ["-o", str(tmp_path / "out"), "--skip-gunc", "--mw1-refinement", "-b", ",".join(dirs)]
    with pytest.raises(SystemExit):
        bin_refinement.main(argv)
    captured = capsys.readouterr()
    said = " ".join((captured.out + captured.err).replace("*", " ").split())
    assert "7 bin sets given, but --mw1-refinement takes at most 6" in said
    assert "The default consolidator (without --mw1-refinement) takes any number" in said
    assert scored_dirs == []


def test_mw1_refinement_works_with_four_sets(tmp_path, scored_dirs, monkeypatch):
    """Every pair and the four-way combination are refined (no set is dropped), and the result
    is laid out like any other run. Binning_refiner runs in process, without its pauses."""
    real_run = bin_refinement.run
    refined = []

    def run(cmd, **kwargs):
        if kwargs.get("tool") != "binning_refiner":
            return real_run(cmd, **kwargs)
        words = shlex.split(cmd)
        refined.append(words[words.index("-o") + 1])
        saved_argv, saved_cwd = sys.argv, os.getcwd()
        monkeypatch.setattr(time, "sleep", lambda _seconds: None)
        sys.argv = ["binning_refiner"] + words[
            words.index("metawrap2.vendor.binning_refiner") + 1 :
        ]
        os.chdir(kwargs["cwd"])
        try:
            runpy.run_module("metawrap2.vendor.binning_refiner", run_name="__main__")
        finally:
            sys.argv = saved_argv
            os.chdir(saved_cwd)
        return 0

    monkeypatch.setattr(bin_refinement, "run", run)
    dirs = _sets(tmp_path, FIVE[:4])
    out = tmp_path / "refined"
    argv = ["-o", str(out), "-t", "1", "-m", "8", "-c", "50", "-x", "10", "--skip-gunc"]
    assert bin_refinement.main(argv + ["--mw1-refinement", "-b", ",".join(dirs)]) == 0
    assert sorted(refined) == sorted(
        ["Refined_AB", "Refined_AC", "Refined_AD", "Refined_BC", "Refined_BD", "Refined_CD"]
        + ["Refined_ABCD"]
    )
    final = os.path.join(str(out), "metawrap_50_10_bins")
    bins = [
        sorted(contig_id(h) for h, _ in iter_fasta(os.path.join(final, f)))
        for f in os.listdir(final)
    ]
    assert bins and all(len(set(b)) == len(b) for b in bins)
