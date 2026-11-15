"""The bin score cache shared between runs (metawrap2.scorecache).

A bin with the same sequences, scored by the same scorer, version, database and options, is not
scored again by a later run. Any change to one of those uses a different cache file. Runs that
share the directory never corrupt it, and a held lock makes a run skip its write instead of
waiting.
"""

from __future__ import annotations

import fcntl
import json
import multiprocessing
import os

import pytest
import test_bin_refinement_golden as golden
import test_golden_end_to_end as e2e
from test_golden_end_to_end import stubs  # noqa: F401 - a fixture

from metawrap2 import checkm, scorecache
from metawrap2.command import runner as _runner
from metawrap2.config import SETTINGS_KEYS, Settings, check_config, load_settings
from metawrap2.constants import BIN_EXTENSION
from metawrap2.modules import _common, bin_refinement

CONFIG = {
    "schema": scorecache.SCHEMA,
    "scorer": "checkm2",
    "tool": "checkm2=1.1.0=x",
    "database": "d",
}


def _write_bin(path, contigs):
    """*contigs* maps contig name to sequence."""
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % item for item in contigs.items())


def _fake_checkm(monkeypatch, scored):
    """Stub run_checkm: scores derived from the bin's size, and a record of what it was given."""

    def fake(bins_dir, **kwargs):
        from metawrap2.io.seqio import fasta_files

        rows = []
        for path in fasta_files(bins_dir):
            name = os.path.splitext(os.path.basename(path))[0]
            scored.append(name)
            total, n50, gc = checkm.bin_shape(path)
            row = [name, "%.2f" % (50 + total % 47), "1.5", "%.3f" % gc, "checkm2", str(n50)]
            rows.append(row + [str(total)] + ([kwargs["binner"]] if kwargs.get("binner") else []))
        out = bins_dir + checkm.STATS_SUFFIX
        checkm.write_stats(rows, out, binner=kwargs.get("binner"))
        return out

    monkeypatch.setattr(checkm, "run_checkm", fake)
    return scored


def _rows(stats):
    with open(stats) as fh:
        return fh.read()


# --- where the cache is and what identifies a bin --------------------------------------------


def test_the_directory_defaults_under_the_user_dir_and_an_empty_value_turns_it_off(monkeypatch):
    monkeypatch.setattr(scorecache, "DEFAULT_DIR", os.path.join("~", ".metawrap2", "score_cache"))
    assert scorecache.cache_dir(None) == os.path.expanduser("~/.metawrap2/score_cache")
    assert scorecache.cache_dir("") is None
    assert scorecache.cache_dir("  ") is None
    assert scorecache.cache_dir("~/elsewhere") == os.path.expanduser("~/elsewhere")


def test_the_config_key_is_known_and_loaded(tmp_path, monkeypatch):
    assert "score_cache_dir" in SETTINGS_KEYS
    assert check_config({"settings": {"score_cache_dir": ""}}) == ([], [])
    assert Settings().score_cache_dir is None
    path = tmp_path / "c.toml"
    path.write_text('[settings]\nscore_cache_dir = "/somewhere/scores"\n')
    assert load_settings(str(path)).score_cache_dir == "/somewhere/scores"
    path.write_text('[settings]\nscore_cache_dir = ""\n')
    assert load_settings(str(path)).score_cache_dir == ""
    path.write_text("[settings]\nthreads = 2\n")
    assert load_settings(str(path)).score_cache_dir is None


def test_a_bin_is_identified_by_its_sequences_not_by_contig_names_or_order(tmp_path):
    _write_bin(tmp_path / "a.fasta", {"c1": "ACGTAC", "c2": "GGGCCC"})
    _write_bin(tmp_path / "b.fasta", {"x9": "GGGCCC", "x1": "ACGTAC"})
    _write_bin(tmp_path / "c.fasta", {"c1": "ACGTAA", "c2": "GGGCCC"})
    a, b, c = (
        scorecache.content_identity(str(tmp_path / n)) for n in ("a.fasta", "b.fasta", "c.fasta")
    )
    assert a == b
    assert a != c
    # The per-run key also hashes each sequence, so a and c differ even though their names match.
    assert checkm.bin_identity(str(tmp_path / "a.fasta")) != checkm.bin_identity(
        str(tmp_path / "c.fasta")
    )


def test_a_database_file_or_directory_is_fingerprinted_and_a_change_shows(tmp_path):
    db = tmp_path / "uniref.dmnd"
    db.write_text("x" * 10)
    first = scorecache.database_identity("CHECKM2DB", str(db))
    assert first and first == scorecache.database_identity("CHECKM2DB", str(db))
    db.write_text("x" * 11)
    assert scorecache.database_identity("CHECKM2DB", str(db)) != first

    data = tmp_path / "checkm_data"
    data.mkdir()
    (data / "taxon_marker_sets.tsv").write_text("t")
    (data / "hmms").mkdir()
    assert scorecache.database_identity("CHECKM_DB", str(data))
    assert scorecache.database_identity("CHECKM_DB", str(tmp_path / "missing")) is None
    assert scorecache.database_identity("CHECKM_DB", "") is None


def _fake_env(root, package="checkm2", version="1.1.0", build="pyh7cba7a3_0", tool="checkm2"):
    (root / "bin").mkdir(parents=True)
    exe = root / "bin" / tool
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    (root / "conda-meta").mkdir()
    (root / "conda-meta" / ("%s-%s-%s.json" % (package, version, build))).write_text("{}")
    (root / "conda-meta" / "python-3.11.0-h0_0.json").write_text("{}")
    return root


def test_the_tool_version_is_read_from_conda_meta_without_running_the_tool(tmp_path, monkeypatch):
    root = _fake_env(tmp_path / "env")
    monkeypatch.setenv("PATH", str(root / "bin"))
    assert scorecache.tool_version("checkm2", "checkm2", None) == "checkm2=1.1.0=pyh7cba7a3_0"

    import metawrap2.config as config_mod

    monkeypatch.setattr(config_mod, "env_prefix", lambda env: str(root) if env == "e" else None)
    assert scorecache.tool_version("checkm2", "checkm2", "e") == "checkm2=1.1.0=pyh7cba7a3_0"
    assert scorecache.tool_version("checkm2", "checkm2", "other") is None
    # checkm-genome is not mistaken for another package that starts the same way.
    assert scorecache.tool_version("checkm-genome", "checkm", "e") is None


def test_the_configuration_covers_scorer_version_database_and_ttable(tmp_path, monkeypatch):
    db = tmp_path / "db.dmnd"
    db.write_text("db")
    monkeypatch.setattr(scorecache, "tool_version", lambda *a: "checkm2=1.1.0=b")
    dbs = {"CHECKM2DB": str(db), "CHECKM_DB": str(db)}

    def path_for(scorer, ttable=None, databases=dbs):
        config, reason = scorecache.scorer_configuration(databases, scorer, ttable, None)
        assert config is not None, reason
        return scorecache.SharedScoreCache(str(tmp_path / "cache"), config).path

    paths = {
        path_for("checkm2"),
        path_for("checkm2", 11),
        path_for("checkm2", 4),
        path_for("checkm1"),
        path_for("checkm1_reduced"),
    }
    assert len(paths) == 5, "each scorer and translation table gets its own file"
    assert path_for("checkm2") == path_for("checkm2"), "the same configuration finds the same file"

    monkeypatch.setattr(scorecache, "tool_version", lambda *a: "checkm2=1.2.0=b")
    assert path_for("checkm2") not in paths, "a new CheckM2 version does not reuse old scores"


@pytest.mark.parametrize(
    "databases, version, expected",
    [
        ({}, "v", "CHECKM2DB is not set"),
        ({"CHECKM2DB": "/no/such/db.dmnd"}, "v", "could not be read"),
        ("db", None, "version could not be read"),
    ],
)
def test_an_unknown_part_of_the_configuration_turns_the_cache_off(
    tmp_path, monkeypatch, databases, version, expected
):
    if databases == "db":
        (tmp_path / "db.dmnd").write_text("db")
        databases = {"CHECKM2DB": str(tmp_path / "db.dmnd")}
    monkeypatch.setattr(scorecache, "tool_version", lambda *a: version)
    config, reason = scorecache.scorer_configuration(databases, "checkm2", None, None)
    assert config is None and expected in reason

    said = []
    shared = scorecache.open_for_run(
        str(tmp_path / "cache"), databases, enabled=True, scorer="checkm2", say=said.append
    )
    assert shared is None
    assert len(said) == 1 and expected in said[0] and "this run only" in said[0]
    assert not (tmp_path / "cache").exists()


def test_the_flag_the_setting_and_a_dry_run_turn_it_off(tmp_path, monkeypatch):
    (tmp_path / "db.dmnd").write_text("db")
    dbs = {"CHECKM2DB": str(tmp_path / "db.dmnd")}
    monkeypatch.setattr(scorecache, "tool_version", lambda *a: "checkm2=1=b")
    said = []

    def opened(setting, enabled=True):
        return scorecache.open_for_run(
            setting, dbs, enabled=enabled, scorer="checkm2", say=said.append
        )

    assert opened(str(tmp_path / "c"), enabled=False) is None
    assert "--no-score-cache" in said[-1]
    assert opened("") is None
    assert "score_cache_dir is empty" in said[-1]
    monkeypatch.setattr(_runner, "dry_run", True)
    assert opened(str(tmp_path / "c")) is None
    monkeypatch.setattr(_runner, "dry_run", False)
    shared = opened(str(tmp_path / "c"))
    assert shared is not None and shared.path in said[-1]


# --- the store -----------------------------------------------------------------------------


def test_scores_written_by_one_run_are_read_by_the_next(tmp_path):
    first = scorecache.SharedScoreCache(str(tmp_path), CONFIG)
    first.put("k1", ["90.0", "1.0", "0.50", "checkm2", "1000", "2000000"])
    assert first.flush()
    assert first.as_record()["added"] == 1

    second = scorecache.SharedScoreCache(str(tmp_path), CONFIG)
    assert second.loaded == 1
    assert second.get("k1") == ["90.0", "1.0", "0.50", "checkm2", "1000", "2000000"]
    assert second.get("k2") is None
    assert (second.hits, second.misses) == (1, 1)

    other = scorecache.SharedScoreCache(str(tmp_path), dict(CONFIG, tool="checkm2=9=x"))
    assert other.loaded == 0 and other.path != first.path


def test_two_writers_merge_instead_of_replacing_each_other(tmp_path):
    a = scorecache.SharedScoreCache(str(tmp_path), CONFIG)
    b = scorecache.SharedScoreCache(str(tmp_path), CONFIG)
    a.put("ka", ["1", "2", "3", "4", "5", "6"])
    b.put("kb", ["1", "2", "3", "4", "5", "7"])
    assert a.flush() and b.flush()
    assert set(scorecache.SharedScoreCache(str(tmp_path), CONFIG).rows) == {"ka", "kb"}


def _writer(directory, index):
    shared = scorecache.SharedScoreCache(directory, CONFIG, lock_timeout=30)
    for j in range(20):
        shared.put("p%d.%d" % (index, j), [str(index), str(j), "0", "x", "1", "2"])
        shared.flush()


def test_concurrent_processes_leave_a_complete_valid_file(tmp_path):
    procs = [
        multiprocessing.get_context("spawn").Process(target=_writer, args=(str(tmp_path), i))
        for i in range(4)
    ]
    for proc in procs:
        proc.start()
    for proc in procs:
        proc.join(120)
        assert proc.exitcode == 0
    shared = scorecache.SharedScoreCache(str(tmp_path), CONFIG)
    assert len(shared.rows) == 80
    assert not [n for n in os.listdir(tmp_path) if n.endswith(".tmp")]


def test_a_held_lock_skips_the_write_without_waiting_long(tmp_path):
    said = []
    shared = scorecache.SharedScoreCache(str(tmp_path), CONFIG, say=said.append, lock_timeout=0.2)
    shared.put("k", ["1", "2", "3", "4", "5", "6"])
    with open(shared.lock_path, "a") as held:
        fcntl.flock(held, fcntl.LOCK_EX)
        assert shared.flush() is False
    assert "locked by another run" in said[-1]
    assert not os.path.exists(shared.path)
    # Still usable by this run, and written at the next flush.
    assert shared.get("k") == ["1", "2", "3", "4", "5", "6"]
    assert shared.flush() is True
    assert scorecache.SharedScoreCache(str(tmp_path), CONFIG).loaded == 1
    assert shared.as_record()["writes_skipped_on_lock_timeout"] == 1


def test_an_unreadable_cache_file_is_ignored_and_replaced(tmp_path):
    shared = scorecache.SharedScoreCache(str(tmp_path), CONFIG)
    with open(shared.path, "w") as fh:
        fh.write('{"config": ')
    again = scorecache.SharedScoreCache(str(tmp_path), CONFIG)
    assert again.loaded == 0
    again.put("k", ["1", "2", "3", "4", "5", "6"])
    assert again.flush()
    with open(again.path) as fh:
        assert json.load(fh)["rows"] == {"k": ["1", "2", "3", "4", "5", "6"]}


def test_a_write_failure_is_said_once_and_does_not_raise(tmp_path):
    blocker = tmp_path / "file"
    blocker.write_text("not a directory")
    said = []
    shared = scorecache.SharedScoreCache(str(blocker / "cache"), CONFIG, say=said.append)
    for key in ("a", "b"):
        shared.put(key, ["1", "2", "3", "4", "5", "6"])
        assert shared.flush() is False
    assert len(said) == 1 and "Could not write" in said[0]


# --- scoring through BinScoreCache ----------------------------------------------------------


def _bins(root, contigs_by_bin):
    root.mkdir()
    for name, contigs in contigs_by_bin.items():
        _write_bin(root / (name + ".fasta"), contigs)
    return str(root)


BIG = {"c%d" % i: ("ACGGT" * (60_000 + 9_000 * i)) for i in range(1, 5)}


def test_a_later_run_reuses_scores_and_writes_the_same_stats(tmp_path, monkeypatch):
    scored = _fake_checkm(monkeypatch, [])
    layout = {"bin_1": {k: BIG[k] for k in ("c1", "c2")}, "bin_2": {"c3": BIG["c3"]}}

    first = _bins(tmp_path / "run1", layout)
    shared = scorecache.SharedScoreCache(str(tmp_path / "cache"), CONFIG)
    cache = checkm.BinScoreCache(str(tmp_path / "run1.json"), shared=shared)
    checkm.run_checkm_cached(first, cache=cache, threads=1, mem_gb=8, env=None)
    assert sorted(scored) == ["bin_1", "bin_2"]
    assert shared.as_record()["added"] == 2

    # A new run in another directory, with renamed contigs and the bins in another order.
    renamed = {
        "bin_A": {"z3": BIG["c3"]},
        "bin_B": {"z2": BIG["c2"], "z1": BIG["c1"]},
        "bin_C": {"z4": BIG["c4"]},
    }
    second = _bins(tmp_path / "run2", renamed)
    shared2 = scorecache.SharedScoreCache(str(tmp_path / "cache"), CONFIG)
    cache2 = checkm.BinScoreCache(shared=shared2)
    del scored[:]
    checkm.run_checkm_cached(second, cache=cache2, threads=1, mem_gb=8, env=None)
    assert scored == ["bin_C"], "only the bin no earlier run scored is sent to CheckM"
    assert (shared2.hits, shared2.misses) == (2, 1)
    assert "2 bin(s) reused from earlier runs" in cache2.summary()

    # The same as scoring the second run's bins from scratch.
    del scored[:]
    cold = _bins(tmp_path / "cold", renamed)
    checkm.run_checkm_cached(cold, cache=checkm.BinScoreCache(), threads=1, mem_gb=8, env=None)
    assert _rows(second + checkm.STATS_SUFFIX) == _rows(cold + checkm.STATS_SUFFIX)


def test_the_same_contig_names_with_other_sequences_are_scored_again(tmp_path, monkeypatch):
    scored = _fake_checkm(monkeypatch, [])
    shared = scorecache.SharedScoreCache(str(tmp_path / "cache"), CONFIG)
    first = _bins(tmp_path / "a", {"bin_1": {"c1": BIG["c1"]}})
    checkm.run_checkm_cached(
        first, cache=checkm.BinScoreCache(shared=shared), threads=1, mem_gb=8, env=None
    )
    other = _bins(tmp_path / "b", {"bin_1": {"c1": BIG["c2"]}})
    shared2 = scorecache.SharedScoreCache(str(tmp_path / "cache"), CONFIG)
    checkm.run_checkm_cached(
        other, cache=checkm.BinScoreCache(shared=shared2), threads=1, mem_gb=8, env=None
    )
    assert scored == ["bin_1", "bin_1"]


def test_a_reused_score_gets_this_runs_binner_column(tmp_path, monkeypatch):
    _fake_checkm(monkeypatch, [])
    shared = scorecache.SharedScoreCache(str(tmp_path / "cache"), CONFIG)
    first = _bins(tmp_path / "a", {"bin_1": {"c1": BIG["c1"]}})
    checkm.run_checkm_cached(
        first, cache=checkm.BinScoreCache(shared=shared), threads=1, mem_gb=8, env=None
    )
    second = _bins(tmp_path / "b", {"bin_9": {"q": BIG["c1"]}, "bin_2": {"c2": BIG["c2"]}})
    shared2 = scorecache.SharedScoreCache(str(tmp_path / "cache"), CONFIG)
    checkm.run_checkm_cached(
        second,
        cache=checkm.BinScoreCache(shared=shared2),
        threads=1,
        mem_gb=8,
        env=None,
        binner="metabat2",
    )
    lines = _rows(second + checkm.STATS_SUFFIX).splitlines()
    assert lines[0].split("\t")[-1] == "binner"
    assert all(len(line.split("\t")) == 8 and line.endswith("\tmetabat2") for line in lines[1:])
    assert any(line.startswith("bin_9\t") for line in lines[1:])


def test_the_plan_leaves_out_bins_an_earlier_run_scored(tmp_path, monkeypatch):
    _fake_checkm(monkeypatch, [])
    shared = scorecache.SharedScoreCache(str(tmp_path / "cache"), CONFIG)
    first = _bins(tmp_path / "a", {"bin_1": {"c1": BIG["c1"]}})
    checkm.run_checkm_cached(
        first, cache=checkm.BinScoreCache(shared=shared), threads=1, mem_gb=8, env=None
    )

    set_a = _bins(tmp_path / "binsA", {"bin_1": {"x": BIG["c1"]}, "bin_2": {"y": BIG["c2"]}})
    set_b = _bins(tmp_path / "binsB", {"bin_1": {"y": BIG["c2"]}, "bin_2": {"w": BIG["c3"]}})
    shared2 = scorecache.SharedScoreCache(str(tmp_path / "cache"), CONFIG)
    plan = checkm.plan_scoring([set_a, set_b], checkm.BinScoreCache(shared=shared2))
    assert [os.path.basename(p) for p in plan[set_a]] == ["bin_2.fasta"]
    assert [os.path.basename(p) for p in plan[set_b]] == ["bin_2.fasta"]


def test_scores_from_the_reduced_tree_fallback_are_not_shared(tmp_path, monkeypatch):
    _fake_checkm(monkeypatch, [])
    shared = scorecache.SharedScoreCache(str(tmp_path / "cache"), dict(CONFIG, scorer="checkm1"))
    first = _bins(tmp_path / "a", {"bin_1": {"c1": BIG["c1"]}})
    monkeypatch.setattr(checkm, "_REDUCED_FALLBACK", {os.path.abspath(first)})
    cache = checkm.BinScoreCache(shared=shared)
    checkm.run_checkm_cached(first, cache=cache, threads=1, mem_gb=8, env=None)
    assert len(cache.rows) == 1, "this run still reuses it"
    assert shared.as_record()["added"] == 0 and not shared.pending


def test_without_a_shared_cache_the_summary_is_unchanged(tmp_path, monkeypatch):
    _fake_checkm(monkeypatch, [])
    cache = checkm.BinScoreCache()
    checkm.run_checkm_cached(
        _bins(tmp_path / "a", {"bin_1": {"c1": BIG["c1"]}}),
        cache=cache,
        threads=1,
        mem_gb=8,
        env=None,
    )
    assert cache.summary() == (
        "1 bin(s): 1 scored, 0 reused from an identical bin, 0 skipped as implausible"
    )


# --- modules ---------------------------------------------------------------------------------


@pytest.fixture
def shared_config(tmp_path, monkeypatch):
    """A config with a CheckM2 database and the shared cache in tmp_path/cache."""
    db = tmp_path / "uniref100.KO.1.dmnd"
    db.write_text("db")
    path = tmp_path / "metawrap2.toml"
    path.write_text(
        '[settings]\nuse_conda_envs = false\nscore_cache_dir = "%s"\n\n[databases]\nCHECKM2DB = "%s"\n'
        % (tmp_path / "cache", db)
    )
    monkeypatch.setenv("METAWRAP2_CONFIG", str(path))
    monkeypatch.setattr(scorecache, "tool_version", lambda *a: "checkm2=1.1.0=b")
    monkeypatch.setattr(_common, "_require_env_manager", lambda env: None)
    return tmp_path / "cache"


def _refine(tmp_path, monkeypatch, name, extra=()):
    calls = []

    def counting(bins_dir, **kwargs):
        calls.extend(n for n in os.listdir(bins_dir) if n.endswith(BIN_EXTENSION))
        return golden._fake_run_checkm(bins_dir, **kwargs)

    monkeypatch.setattr(checkm, "run_checkm", counting)
    root = tmp_path / name
    root.mkdir()
    dirs = golden._write_bin_sets(root, "ABC")
    out = root / "refined"
    argv = ["-o", str(out), "-t", "1", "-m", "8", "-x", "10", "--skip-gunc"]
    argv += ["-A", dirs["A"], "-B", dirs["B"], "-C", dirs["C"]] + list(extra)
    assert bin_refinement.main(argv) == 0
    with open(out / "run_config.json") as fh:
        record = json.load(fh)["score_cache"]
    return str(out), calls, record


def test_a_second_refinement_with_other_thresholds_scores_nothing(
    tmp_path, monkeypatch, shared_config
):
    _out, calls, record = _refine(tmp_path, monkeypatch, "first", ["-c", "50"])
    assert calls and record["status"] == "used" and record["added"] > 0
    assert record["configuration"]["scorer"] == "checkm2"
    assert record["configuration"]["ttable"] == "both"

    out, calls, record = _refine(tmp_path, monkeypatch, "second", ["-c", "60"])
    assert calls == [], "every candidate was scored by the first run"
    assert record["hits"] > 0 and record["misses"] == 0 and record["added"] == 0
    with open(os.path.join(out, "metawrap2.log")) as fh:
        assert "reused from earlier runs" in fh.read()

    # The same result as a run that scores everything itself.
    cold, calls, record = _refine(tmp_path, monkeypatch, "cold", ["-c", "60", "--no-score-cache"])
    assert calls and record == {"status": "off", "reason": "--no-score-cache"}

    # Leaving out the stats of the temporary sets CheckM was run on, which only a run that
    # scores leaves behind.
    def kept(root):
        return {k: v for k, v in golden._snapshot(root).items() if ".to_score." not in k}

    assert kept(out) == kept(cold)


def test_another_translation_table_does_not_reuse_scores(tmp_path, monkeypatch):
    """The translation-table choice is part of the configuration, so it selects another file."""
    db = tmp_path / "uniref.dmnd"
    db.write_text("x")
    monkeypatch.setattr(scorecache, "tool_version", lambda *a, **k: "checkm2=1.0.2=py_0")
    both, _ = scorecache.scorer_configuration({"CHECKM2DB": str(db)}, "checkm2", None, None)
    t11, _ = scorecache.scorer_configuration({"CHECKM2DB": str(db)}, "checkm2", 11, None)
    assert both["ttable"] == "both" and t11["ttable"] == "11"
    assert both != t11


def test_every_scoring_module_has_the_flag():
    from metawrap2.commands import score_cmd
    from metawrap2.modules import binning, reassemble_bins

    assert (
        bin_refinement._parse_args(["-o", "o", "-A", "a", "--no-score-cache"]).score_cache is False
    )
    assert bin_refinement._parse_args(["-o", "o", "-A", "a"]).score_cache is True
    assert score_cmd._parse_args(["-o", "o", "d", "--no-score-cache"]).score_cache is False
    assert (
        reassemble_bins._parse_args(
            ["-o", "o", "-b", "b", "-1", "r1", "-2", "r2", "--no-score-cache"]
        ).score_cache
        is False
    )
    assert (
        binning._parse_args(
            ["-o", "o", "-a", "a.fa", "--run-checkm", "--no-score-cache", "r_1.fastq", "r_2.fastq"]
        ).score_cache
        is False
    )


def test_bin_refinement_and_scorecache_name_scorers_the_same():
    from types import SimpleNamespace

    for checkm2 in (True, False):
        for full_tree in (True, False):
            args = SimpleNamespace(checkm2=checkm2, full_tree=full_tree)
            assert bin_refinement.scorer_name(args) == scorecache.scorer_name(checkm2, full_tree)


def test_the_benchmark_measures_arms_without_the_shared_cache():
    from metawrap2.benchmark import runner

    assert runner._scores_bins("bin_refinement", [])
    assert not runner._scores_bins("bin_refinement", ["--skip-checkm"])
    assert runner._scores_bins("reassemble_bins", [])
    assert runner._scores_bins("binning", ["--run-checkm"])
    assert not runner._scores_bins("binning", [])
    assert not runner._scores_bins("assembly", [])


def _scoring_config(tmp_path, gunc=False):
    """A config file with a CheckM2 database (and GUNC's), conda envs off and the cache on."""
    db = tmp_path / "uniref100.KO.1.dmnd"
    db.write_text("db")
    text = '[settings]\nuse_conda_envs = false\nscore_cache_dir = "%s"\n\n[databases]\n' % (
        tmp_path / "cache"
    )
    text += 'CHECKM2DB = "%s"\n' % db
    if gunc:
        (tmp_path / "gunc.dmnd").write_text("g")
        text += 'GUNC_DB = "%s"\n' % (tmp_path / "gunc.dmnd")
    path = tmp_path / "scoring.toml"
    path.write_text(text)
    return str(path)


@pytest.mark.usefixtures("stubs")
def test_binning_run_checkm_reuses_scores_from_an_earlier_run(tmp_path, monkeypatch):
    from metawrap2.modules import binning

    monkeypatch.setattr(scorecache, "tool_version", lambda *a: "checkm2=1.1.0=b")
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    config = _scoring_config(tmp_path)
    scored = _fake_checkm(monkeypatch, [])
    asm = tmp_path / "asm.fa"
    asm.write_text(e2e.ASSEMBLY)
    for suffix in ("_1", "_2"):
        (tmp_path / ("s%s.fastq" % suffix)).write_text(
            "".join("@r%d/1\nACGT\n+\nIIII\n" % i for i in e2e.CONTIGS)
        )

    def run(name, *extra):
        out = tmp_path / name
        argv = ["-a", str(asm), "-o", str(out), "--metabat2", "--maxbin2", "-t", "1"]
        argv += ["--run-checkm", "--config", config] + list(extra)
        argv += [str(tmp_path / "s_1.fastq"), str(tmp_path / "s_2.fastq")]
        assert binning.main(argv) == 0
        with open(out / "run_config.json") as fh:
            return out, json.load(fh)["score_cache"]

    first, record = run("first")
    # The stub contigs share one sequence, so some bins have the same contents.
    assert scored and record["added"] >= 1
    count = len(scored)
    second, record = run("second")
    assert len(scored) == count, "no bin was scored again"
    assert record["hits"] == 5 and record["misses"] == 0  # 2 metaBAT2 bins, 3 MaxBin2 bins
    for name in ("metabat2_bins", "maxbin2_bins"):
        stats = name + checkm.STATS_SUFFIX
        assert _rows(str(first / stats)) == _rows(str(second / stats))

    _third, record = run("third", "--no-score-cache")
    assert len(scored) == count + 5 and record["status"] == "off"  # every bin scored


def test_reassemble_bins_scores_the_kept_versions_only_once(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from metawrap2.modules import reassemble_bins

    scored = _fake_checkm(monkeypatch, [])
    shared = scorecache.SharedScoreCache(str(tmp_path / "cache"), CONFIG)
    bins = _bins(
        tmp_path / "reassembled_bins",
        {"bin_1.orig": {"NODE_1": BIG["c1"]}, "bin_1.strict": {"NODE_1": BIG["c2"]}},
    )
    args = SimpleNamespace(threads=1, memory=8)
    first = _rows(reassemble_bins._score(args, bins, "env", shared))
    assert sorted(scored) == ["bin_1.orig", "bin_1.strict"], "same contig name, other sequence"
    os.remove(os.path.join(bins, "bin_1.orig.fasta"))
    second = _rows(reassemble_bins._score(args, bins, "env", shared))
    assert len(scored) == 2, "the kept version was not scored again"
    assert second.splitlines()[1] in first.splitlines()

    # Without the shared cache, CheckM2 runs as before.
    reassemble_bins._score(args, bins, "env", None)
    assert len(scored) == 3


def test_score_reuses_scores_from_an_earlier_score_run(tmp_path, monkeypatch):
    from metawrap2.commands import score_cmd

    monkeypatch.setattr(scorecache, "tool_version", lambda *a: "checkm2=1.1.0=b")
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    config = _scoring_config(tmp_path)
    scored = _fake_checkm(monkeypatch, [])
    bins = _bins(tmp_path / "mags", {"bin_1": {"c1": BIG["c1"]}, "bin_2": {"c2": BIG["c2"]}})

    def run(name):
        out = tmp_path / name
        argv = ["-o", str(out), "-t", "1", "-m", "8", "--no-gunc", "--config", config, bins]
        assert score_cmd.main(argv) == 0
        with open(out / "run_config.json") as fh:
            return json.load(fh)["score_cache"]

    assert run("first")["added"] == 2
    assert run("second")["hits"] == 2
    assert len(scored) == 2


def test_score_into_an_existing_run_says_how_to_continue_and_exits_1(tmp_path, monkeypatch, capsys):
    from metawrap2.commands import score_cmd

    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    out = tmp_path / "scored"
    out.mkdir()
    (out / "run_config.json").write_text('{"command_line": ["metawrap2", "score"]}')
    bins = _bins(tmp_path / "mags", {"bin_1": {"c1": BIG["c1"]}})
    with pytest.raises(SystemExit) as stopped:
        score_cmd.main(["-o", str(out), "--no-gunc", bins])
    assert stopped.value.code == 1
    from _banner import unbanner

    said = " ".join(unbanner(capsys.readouterr().out).split())
    assert str(out) in said
    assert "run the same command again with --resume" in said
    assert "checkm_scores.checkm2.json" in said
    assert "add --force" in said and "choose a new -o" in said


def test_config_show_names_the_score_cache(tmp_path, capsys):
    from metawrap2.commands import config_cmd

    path = tmp_path / "c.toml"
    path.write_text('[settings]\nscore_cache_dir = "%s"\n' % (tmp_path / "scores"))
    assert config_cmd.main(["show", "--config", str(path)]) == 0
    assert "score cache:    %s" % (tmp_path / "scores") in capsys.readouterr().out
    path.write_text('[settings]\nscore_cache_dir = ""\n')
    assert config_cmd.main(["show", "--config", str(path)]) == 0
    assert "score cache:    <off>" in capsys.readouterr().out


def test_bins_with_the_same_contig_names_but_different_sequences_get_different_keys(tmp_path):
    """Two single-contig genomes both named `>1` must not share one cached score."""
    from metawrap2 import checkm

    a = tmp_path / "a.fasta"
    b = tmp_path / "b.fasta"
    a.write_text(">1\nACGTACGT\n")
    b.write_text(">1\nTTTTGGGG\n")
    assert checkm.bin_identity(str(a)) != checkm.bin_identity(str(b))
    c = tmp_path / "c.fasta"
    c.write_text(">2\nTTTT\n>1\nACGTACGT\n")
    d = tmp_path / "d.fasta"
    d.write_text(">1\nACGTACGT\n>2\nTTTT\n")
    assert checkm.bin_identity(str(c)) == checkm.bin_identity(str(d))
