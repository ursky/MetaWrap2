"""Tests for bin_refinement's default consolidator (metawrap2) and its algorithm in metawrap2.consolidation.

Also covers how --consolidator, --mw1-refinement and the options only metaWRAP 1's method uses
select the method.

The algorithm tests use small synthetic bins and made-up scores with a known right answer. The
module tests run bin_refinement with CheckM replaced by a stand-in that scores a bin against the
genome it mostly comes from, so a union that joins two halves of a genome scores as complete.
"""

from __future__ import annotations

import itertools
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest

from metawrap2 import checkm as checkm_mod
from metawrap2 import command, consolidation, scorecache
from metawrap2.constants import BIN_EXTENSION, CONTAMINATION_WEIGHT, STATS_SUFFIX
from metawrap2.io.seqio import contig_id, iter_fasta
from metawrap2.modules import _common, bin_refinement

#: Every contig is 100 kb, so a two-contig candidate clears the 200 kb floor and nothing reaches the
#: 20 Mb ceiling. The length window only matters in the test about it.
LEN = 100_000


def lengths_for(bins):
    return {c: LEN for contigs in bins.values() for c in contigs}


def as_bins(spec):
    return {name: frozenset(contigs) for name, contigs in spec.items()}


def pool_of(bins, operations=consolidation.OPERATIONS, lengths=None, **kwargs):
    return consolidation.build_candidates(bins, lengths or lengths_for(bins), operations, **kwargs)


def scored(pool, score_fn):
    """``{signature: (completeness, contamination)}`` from a function of kind and parents."""
    return {key: score_fn(c.kind, c.parents) for key, c in pool.items()}


def total(selection):
    return sum(s.score for s in selection.chosen)


def contig_sets(pool, kind=None):
    return {c.contigs for c in pool.values() if kind is None or c.kind == kind}


#: A chain: a overlaps b, b overlaps c, and a and c share nothing.
CHAIN = as_bins(
    {"a/bin1": {"c1", "c2", "c3"}, "b/bin1": {"c3", "c4", "c5"}, "c/bin1": {"c5", "c6", "c7"}}
)


def chain_scores(kind, parents):
    """The ends score 80 each, the middle 70 and any derived candidate 85 (weight 5, no
    contamination)."""
    if kind == consolidation.ORIGINAL:
        return (70.0, 0.0) if parents[0] == "b/bin1" else (80.0, 0.0)
    return (85.0, 0.0)


# --- candidate generation -----------------------------------------------------------------


def test_a_union_of_two_halves_of_a_genome_is_a_candidate():
    bins = as_bins({"a/bin1": {"c1", "c2", "c3"}, "b/bin1": {"c3", "c4", "c5"}})
    assert frozenset({"c1", "c2", "c3", "c4", "c5"}) in contig_sets(pool_of(bins), "union")


def test_an_intersection_only_pool_has_no_candidate_larger_than_an_input_bin():
    bins = as_bins({"a/bin1": {"c1", "c2", "c3", "c4"}, "b/bin1": {"c3", "c4", "c5", "c6"}})
    pool = pool_of(bins, ("intersection",))
    assert all(any(c <= b for b in bins.values()) for c in contig_sets(pool))
    assert contig_sets(pool, "intersection") == {frozenset({"c3", "c4"})}


def test_differences_are_candidates():
    bins = as_bins({"a/bin1": {"c1", "c2", "c3", "c4"}, "b/bin1": {"c3", "c4"}})
    pool = pool_of(bins, ("difference",))
    assert contig_sets(pool, "difference") == {frozenset({"c1", "c2"})}
    (diff,) = [c for c in pool.values() if c.kind == "difference"]
    # The bin the others are subtracted from comes first, then the clique.
    assert diff.parents == ("a/bin1", "a/bin1", "b/bin1")


def test_only_cliques_are_combined_so_a_chain_gives_no_union_of_its_ends():
    pool = pool_of(CHAIN, ("union",))
    unions = sorted(c.parents for c in pool.values() if c.kind == "union")
    assert unions == [("a/bin1", "b/bin1"), ("b/bin1", "c/bin1")]
    assert frozenset({"c1", "c2", "c3", "c5", "c6", "c7"}) not in contig_sets(pool)


def test_a_three_way_clique_is_combined():
    bins = as_bins(
        {"a/1": {"c1", "c2", "c3", "c6"}, "b/1": {"c2", "c3", "c7"}, "c/1": {"c1", "c3", "c8"}}
    )
    pool = pool_of(bins, ("union",))
    assert ("a/1", "b/1", "c/1") in [c.parents for c in pool.values()]


def test_cliques_match_filtering_every_combination():
    rng = random.Random(3)
    names = ["n%02d" % i for i in range(14)]
    adjacency = {n: set() for n in names}
    for a, b in itertools.combinations(names, 2):
        if rng.random() < 0.4:
            adjacency[a].add(b)
            adjacency[b].add(a)
    for size in (2, 3, 4):
        expected = [
            combo
            for combo in itertools.combinations(names, size)
            if all(y in adjacency[x] for x, y in itertools.combinations(combo, 2))
        ]
        assert list(consolidation.cliques(names, adjacency, size)) == expected


def test_the_length_window_drops_derived_candidates_and_keeps_input_bins():
    bins = as_bins({"a/bin1": {"c1", "c2"}, "b/bin1": {"c2", "c3"}})
    lengths = {"c1": 50_000, "c2": 50_000, "c3": 50_000}
    pool = pool_of(bins, lengths=lengths)
    assert sorted(c.kind for c in pool.values()) == [consolidation.ORIGINAL, consolidation.ORIGINAL]

    huge = {"c1": 8_000_000, "c2": 8_000_000, "c3": 8_000_000}
    pool = pool_of(bins, ("union",), lengths=huge)
    assert "union" not in {c.kind for c in pool.values()}, "24 Mb is above the 20 Mb ceiling"
    assert consolidation.MIN_LENGTH == checkm_mod.MIN_PLAUSIBLE_BIN_BP == 200_000
    assert consolidation.MAX_LENGTH == checkm_mod.MAX_PLAUSIBLE_BIN_BP == 20_000_000


def test_bins_sharing_under_a_tenth_of_the_smaller_bin_are_not_combined():
    bins = as_bins(
        {
            "a/bin1": {"c%d" % i for i in range(20)},
            "b/bin1": {"c0"} | {"d%d" % i for i in range(20)},
        }
    )
    lengths = lengths_for(bins)
    assert consolidation.MIN_OVERLAP_FRACTION == 0.10
    assert consolidation.overlap_graph(bins, lengths) == {}
    assert consolidation.overlap_graph(bins, lengths, min_overlap=0.0)


def test_bins_more_than_four_times_apart_in_length_are_not_combined():
    bins = as_bins({"a/big": {"c%d" % i for i in range(100)}, "b/small": {"c0", "c1", "c2"}})
    lengths = lengths_for(bins)
    assert consolidation.MAX_LENGTH_RATIO == 4.0
    assert consolidation.overlap_graph(bins, lengths) == {}
    assert consolidation.overlap_graph(bins, lengths, max_length_ratio=0)


def test_identical_contig_sets_are_kept_once_under_the_first_kind():
    """Two binners' identical bins are one candidate; a derived set equal to an input bin stays
    an input bin, so it gets the input-bin tiebreak and is not checked by the lookahead."""
    bins = as_bins(
        {
            "a/bin1": {"c1", "c2", "c3"},
            "b/bin1": {"c1", "c2", "c3"},
            "b/bin2": {"c1", "c2", "c3", "c4", "c5"},
        }
    )
    pool = pool_of(bins)
    assert len({c.contigs for c in pool.values()}) == len(pool)
    first = [c for c in pool.values() if c.contigs == frozenset({"c1", "c2", "c3"})]
    assert [(c.kind, c.parents) for c in first] == [(consolidation.ORIGINAL, ("a/bin1",))]
    # b/bin2 minus a/bin1 is {c4, c5}; the union with b/bin2 is b/bin2 itself.
    assert frozenset({"c4", "c5"}) in contig_sets(pool, "difference")
    assert "union" not in {c.kind for c in pool.values()}


def test_an_unknown_operation_is_refused():
    with pytest.raises(ValueError, match="unknown operation"):
        pool_of(CHAIN, ("merge",))


# --- selection ----------------------------------------------------------------------------


def test_selected_bins_share_no_contig():
    bins = as_bins(
        {"a/bin1": {"c1", "c2", "c3"}, "b/bin1": {"c2", "c3", "c4"}, "c/bin1": {"c7", "c8"}}
    )
    pool = pool_of(bins)
    result = consolidation.select(pool, {k: (95.0, 1.0) for k in pool}, 50, 10)
    seen = set()
    for picked in result.chosen:
        assert not picked.candidate.contigs & seen
        seen |= picked.candidate.contigs
    assert result.blocked > 0


def test_a_union_scoring_below_its_parent_is_not_selected():
    bins = as_bins({"a/bin1": {"c1", "c2"}, "b/bin1": {"c2", "c3"}})
    pool = pool_of(bins, ("union",))

    def score_fn(kind, parents):
        if kind == consolidation.ORIGINAL:
            return (90.0, 1.0) if parents[0] == "a/bin1" else (60.0, 1.0)
        return (95.0, 4.0)  # 75, below a/bin1's 85

    result = consolidation.select(pool, scored(pool, score_fn), 50, 10)
    assert [s.candidate.parents for s in result.chosen] == [("a/bin1",)]


def test_a_union_scoring_above_both_parents_is_selected():
    bins = as_bins({"a/bin1": {"c1", "c2"}, "b/bin1": {"c2", "c3"}})
    pool = pool_of(bins, ("union",))

    def score_fn(kind, _parents):
        return (50.0, 1.0) if kind == consolidation.ORIGINAL else (99.0, 1.0)

    result = consolidation.select(pool, scored(pool, score_fn), 40, 10)
    assert [s.candidate.kind for s in result.chosen] == ["union"]


def test_without_lookahead_one_union_displaces_two_disjoint_bins():
    """On the chain, the union of a and b holds the contig c needs, so it blocks a and c: 85
    against the 160 the input bins alone reach."""
    inter = pool_of(CHAIN, ("intersection",))
    union = pool_of(CHAIN, ("union",))
    base = consolidation.select(inter, scored(inter, chain_scores), 50, 10, lookahead=False)
    greedy = consolidation.select(union, scored(union, chain_scores), 50, 10, lookahead=False)
    assert sorted(s.candidate.parents for s in base.chosen) == [("a/bin1",), ("c/bin1",)]
    assert total(base) == 160
    assert [s.candidate.kind for s in greedy.chosen] == ["union"]
    assert total(greedy) == 85


def test_the_lookahead_refuses_that_union_and_keeps_the_two_bins():
    union = pool_of(CHAIN, ("union",))
    guarded = consolidation.select(union, scored(union, chain_scores), 50, 10, lookahead=True)
    assert total(guarded) == 160
    assert {s.candidate.kind for s in guarded.chosen} == {consolidation.ORIGINAL}
    assert guarded.rejected_by_lookahead == 2


def test_the_lookahead_packs_its_alternative_over_input_bins_and_intersections_only():
    """With every operation the pool also holds a competing union worth 85. Counting it as the
    alternative would let the first union through."""
    pool = pool_of(CHAIN)
    assert {"union", "difference"} <= {c.kind for c in pool.values()}

    def score_fn(kind, parents):
        if kind == "difference":
            return (10.0, 0.0)  # below -c, not considered
        return chain_scores(kind, parents)

    guarded = consolidation.select(pool, scored(pool, score_fn), 50, 10)
    assert total(guarded) == 160


def test_the_lookahead_takes_a_union_that_outscores_what_it_blocks():
    """The union of a and b blocks a and c. It is refused at 99 against 60 + 60 and taken at 99
    against 45 + 45."""
    union = pool_of(CHAIN, ("union",))

    def scores_with(original):
        def score_fn(kind, parents):
            if kind == consolidation.ORIGINAL:
                return (original, 0.0)
            return (99.0, 0.0) if parents == ("a/bin1", "b/bin1") else (41.0, 0.0)

        return scored(union, score_fn)

    refused = consolidation.select(union, scores_with(60.0), 40, 10)
    assert "union" not in {s.candidate.kind for s in refused.chosen}
    assert total(refused) == 120
    taken = consolidation.select(union, scores_with(45.0), 40, 10)
    assert taken.chosen[0].candidate.parents == ("a/bin1", "b/bin1")
    assert total(taken) == 99


def test_the_contamination_weight_changes_the_winner():
    bins = as_bins({"a/bin1": {"c1", "c2"}, "b/bin1": {"c1", "c2", "c3"}})
    pool = pool_of(bins, ())

    def score_fn(_kind, parents):
        return (80.0, 1.0) if parents[0] == "a/bin1" else (88.0, 2.8)

    scores = scored(pool, score_fn)
    at_5 = consolidation.select(pool, scores, 50, 10, weight=5)
    at_2 = consolidation.select(pool, scores, 50, 10, weight=2)
    assert [s.candidate.parents for s in at_5.chosen] == [("a/bin1",)]  # 75 against 74
    assert [s.candidate.parents for s in at_2.chosen] == [("b/bin1",)]  # 78 against 82.4
    assert consolidation.select(pool, scores, 50, 10).chosen == at_5.chosen
    assert CONTAMINATION_WEIGHT == 5.0


def test_the_thresholds_are_strict_as_in_the_default_consolidator():
    bins = as_bins({"a/bin1": {"c1", "c2"}, "b/bin1": {"c5", "c6"}, "c/bin1": {"c8", "c9"}})
    pool = pool_of(bins, ())
    scores = scored(
        pool,
        lambda _k, p: {"a/bin1": (70.0, 1.0), "b/bin1": (90.0, 10.0), "c/bin1": (70.1, 9.9)}[p[0]],
    )
    result = consolidation.select(pool, scores, 70, 10)
    assert [s.candidate.parents for s in result.chosen] == [("c/bin1",)]


def test_ties_go_to_higher_n50_then_input_bins_then_signature():
    bins = as_bins({"a/bin1": {"c1", "c2"}, "b/bin1": {"c2", "c3"}})
    pool = pool_of(bins, ())
    flat = {k: (90.0, 1.0) for k in pool}

    def winner(lengths):
        result = consolidation.select(pool, flat, 50, 10, lengths=lengths)
        return [s.candidate.parents[0] for s in result.chosen]

    # The two bins overlap and score the same, so N50 decides which one is taken.
    assert winner({"c1": 300_000, "c2": 100_000, "c3": 150_000}) == ["a/bin1"]
    assert winner({"c1": 150_000, "c2": 100_000, "c3": 300_000}) == ["b/bin1"]
    # With no lengths the signature decides.
    assert winner(None) == [pool[min(pool)].parents[0]]

    # An input bin and a union with the same score and N50: the input bin is taken.
    lengths = {"x1": 150_000, "x2": 150_000, "y1": 150_000}
    bins = as_bins({"a/1": {"x1", "x2"}, "b/1": {"x2", "y1"}})
    pool = pool_of(bins, ("union",), lengths=lengths)
    keep = {k: (90.0, 1.0) for k, c in pool.items() if c.parents != ("b/1",)}
    assert {c.kind for k, c in pool.items() if k in keep} == {consolidation.ORIGINAL, "union"}
    result = consolidation.select(pool, keep, 50, 10, lengths=lengths, lookahead=False)
    assert [(s.candidate.kind, s.candidate.parents) for s in result.chosen] == [
        (consolidation.ORIGINAL, ("a/1",))
    ]


def test_input_order_does_not_change_the_pool_or_the_selection():
    rng = random.Random(11)
    contigs = ["k%03d" % i for i in range(60)]
    spec = {}
    for binner in "abc":
        shuffled = contigs[:]
        rng.shuffle(shuffled)
        for i in range(0, 60, 12):
            spec["%s/bin%d" % (binner, i)] = set(shuffled[i : i + 12 + rng.randint(0, 6)])
    bins = as_bins(spec)
    lengths = {c: 50_000 + 1000 * i for i, c in enumerate(contigs)}
    forward = consolidation.build_candidates(bins, lengths, min_overlap=0)
    backward = consolidation.build_candidates(
        dict(reversed(list(bins.items()))), lengths, min_overlap=0
    )
    assert list(forward.items()) == list(backward.items())
    scores = {k: (50 + int(k[:4], 16) % 50, int(k[4:6], 16) % 8) for k in forward}
    a = consolidation.select(forward, scores, 50, 10, lengths=lengths)
    b = consolidation.select(backward, scores, 50, 10, lengths=lengths)
    assert a == b and a.chosen


#: Builds a pool and selects from it, printing the result; run under two PYTHONHASHSEED values.
_DETERMINISM_SCRIPT = """
import json, random
from metawrap2 import consolidation
rng = random.Random(5)
contigs = ["contig_%d" % i for i in range(80)]
bins = {}
for binner in ("metabat2", "semibin2", "comebin"):
    order = contigs[:]
    rng.shuffle(order)
    for i in range(0, 80, 16):
        bins["%s/bin_%d" % (binner, i)] = frozenset(order[i : i + 16 + rng.randint(0, 8)])
lengths = {c: 40000 + 977 * i for i, c in enumerate(contigs)}
pool = consolidation.build_candidates(bins, lengths, min_overlap=0)
scores = {k: (50 + int(k[:4], 16) % 50, int(k[4:6], 16) % 8) for k in pool}
result = consolidation.select(pool, scores, 50, 10, lengths=lengths)
print(json.dumps({
    "pool": [[c.signature, c.kind, list(c.parents)] for c in pool.values()],
    "chosen": [s.candidate.signature for s in result.chosen],
    "counts": [result.blocked, result.rejected_by_lookahead],
}))
"""


def test_the_result_does_not_depend_on_the_hash_seed():
    src = str(Path(consolidation.__file__).parents[1])
    outputs = []
    for seed in ("0", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed, PYTHONPATH=src)
        done = subprocess.run(
            [sys.executable, "-c", _DETERMINISM_SCRIPT],
            capture_output=True,
            text=True,
            check=True,
            env=env,
        )
        outputs.append(json.loads(done.stdout))
    assert outputs[0] == outputs[1]
    assert outputs[0]["chosen"]
    assert {kind for _sig, kind, _p in outputs[0]["pool"]} >= {"original", "union"}


# --- bin_refinement's default consolidator ---------------------------------------------------

CONTIGS = {"c%02d" % i: ("ACGT" * (50000 + i * 2500)) for i in range(1, 13)}
#: Two genomes. Both binners split the first one, at different places; the union of their halves
#: is the whole genome, which the default consolidator cannot produce.
GENOMES = {"g1": ["c01", "c02", "c03", "c04", "c05", "c06"], "g2": ["c07", "c08", "c09", "c10"]}
GENOME_OF = {c: g for g, members in GENOMES.items() for c in members}
BIN_SETS = {
    "A": {"bin.1": ["c01", "c02", "c03", "c04"], "bin.2": ["c07", "c08", "c09", "c10"]},
    "B": {"bin.1": ["c03", "c04", "c05", "c06"], "bin.2": ["c07", "c08", "c09", "c11"]},
}


def _truth_scores(contigs):
    """Completeness against the bin's main genome and contamination from the other contigs."""
    counts = {}
    for c in contigs:
        counts[GENOME_OF.get(c, "other")] = counts.get(GENOME_OF.get(c, "other"), 0) + 1
    main = max(sorted(g for g in counts if g != "other") or ["other"], key=counts.get)
    if main == "other":
        return 0.0, 100.0
    size = len(GENOMES[main])
    return 100.0 * counts[main] / size, 100.0 * (len(contigs) - counts[main]) / size / 4


def _fake_run_checkm(bins_dir, **kwargs):
    rows = []
    for name in sorted(os.listdir(bins_dir)):
        if not name.endswith(BIN_EXTENSION):
            continue
        path = os.path.join(bins_dir, name)
        records = list(iter_fasta(path))
        comp, cont = _truth_scores([contig_id(h) for h, _ in records])
        rows.append(
            [
                os.path.splitext(name)[0],
                "%.2f" % comp,
                "%.2f" % cont,
                "0.500",
                "Bacteria",
                str(max(len(s) for _, s in records)),
                str(sum(len(s) for _, s in records)),
            ]
        )
    stats_path = bins_dir + STATS_SUFFIX
    checkm_mod.write_stats(rows, stats_path)
    return stats_path


@pytest.fixture
def scored_dirs(monkeypatch):
    """Fake CheckM; returns the list of bin directories it was run on."""
    calls = []

    def counting(bins_dir, **kwargs):
        calls.append(os.path.basename(bins_dir))
        return _fake_run_checkm(bins_dir, **kwargs)

    monkeypatch.setattr(bin_refinement._checkm, "run_checkm", counting)
    monkeypatch.setattr(checkm_mod, "run_checkm", counting)
    monkeypatch.setattr(_common, "_require_env_manager", lambda env: None)
    return calls


def _inputs(root):
    dirs = []
    for letter in "AB":
        directory = root / ("binner_%s" % letter)
        directory.mkdir(parents=True, exist_ok=True)
        for bin_name, contigs in BIN_SETS[letter].items():
            with open(directory / ("%s.fa" % bin_name), "w") as fh:
                fh.writelines(">%s\n%s\n" % (c, CONTIGS[c]) for c in contigs)
        dirs.append(str(directory))
    return dirs


def _run(root, extra=(), resume=False, out_name="refined"):
    dirs = _inputs(root)
    out = root / out_name
    argv = ["-o", str(out), "-t", "1", "-m", "8", "-c", "50", "-x", "10", "--skip-gunc"]
    argv += ["-A", dirs[0], "-B", dirs[1], *extra]
    command.set_resume(resume)
    try:
        assert bin_refinement.main(argv) == 0
    finally:
        command.set_resume(False)
    return str(out)


def _final_bins(out):
    final = os.path.join(out, "metawrap_50_10_bins")
    return {
        name: sorted(contig_id(h) for h, _ in iter_fasta(os.path.join(final, name)))
        for name in sorted(os.listdir(final))
    }


def _layout(out):
    """Every path in the output except work_files and .metawrap2."""
    found = set()
    for root, dirs, files in os.walk(out):
        rel = os.path.relpath(root, out)
        if rel.split(os.sep)[0] in ("work_files", ".metawrap2"):
            continue
        found.update(os.path.normpath(os.path.join(rel, f)) for f in files + dirs)
    # Bin numbers differ between consolidators; the names and files around them do not.
    return {re.sub(r"bin_\d+", "bin_N", p) for p in found}


def test_the_default_joins_two_halves_of_a_genome_that_mw1_refinement_cannot(tmp_path, scored_dirs):
    out = _run(tmp_path / "default")
    final = _final_bins(out)
    assert sorted(final.values()) == [GENOMES["g1"], GENOMES["g2"]]

    work = os.path.join(out, "work_files")
    with open(os.path.join(work, bin_refinement.CONSOLIDATED_BINS_FILE)) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    assert rows[0] == [
        "bin",
        "kind",
        "parents",
        "completeness",
        "contamination",
        "score",
        "candidate",
    ]
    assert [(r[0], r[1], r[2], r[3]) for r in rows[1:]] == [
        # Both score 100; the second genome's contigs are longer, so its N50 is higher.
        ("bin_001", "original", "binsA/bin.2", "100.00"),
        ("bin_002", "union", "binsA/bin.1,binsB/bin.1", "100.00"),
    ]
    # Derived candidates are scored as their own bin sets, beside the input sets.
    assert {"binsA", "binsB", "bins_candidates_union"} <= set(scored_dirs)
    assert not any(d.startswith("binsAB") for d in os.listdir(work))
    assert not os.path.isdir(os.path.join(work, "binsM"))

    # The final stats table has the default's header and the same columns.
    with open(os.path.join(out, "metawrap_50_10_bins" + STATS_SUFFIX)) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        body = [line.split("\t") for line in fh]
    assert header[: len(checkm_mod.STATS_HEADER)] == list(checkm_mod.STATS_HEADER)
    assert [r[0] for r in body] == ["bin_001", "bin_002"] and body[0][1] == "100.00"

    # --skip-refinement keeps the metaWRAP 1 run fast; its bins are still subsets of input bins.
    mw1 = _run(tmp_path / "mw1", ["--mw1-refinement", "--skip-refinement"])
    assert max(len(c) for c in _final_bins(mw1).values()) < len(GENOMES["g1"])
    assert _layout(out) == _layout(mw1)
    assert not os.path.isfile(os.path.join(mw1, "work_files", bin_refinement.CONSOLIDATOR_FILE))


def test_a_resumed_run_rescores_nothing_and_gives_the_same_bins(tmp_path, scored_dirs):
    out = _run(tmp_path)
    before = _final_bins(out)
    del scored_dirs[:]
    _run(tmp_path, resume=True)
    assert scored_dirs == []
    assert _final_bins(out) == before


def test_the_weight_reaches_the_selection_and_the_record(tmp_path, scored_dirs):
    out = _run(tmp_path, ["--contamination-weight", "2"])
    with open(os.path.join(out, "metawrap2.log")) as fh:
        assert "at contamination weight 2;" in fh.read()
    with open(os.path.join(out, "run_config.json")) as fh:
        text = fh.read()
    assert '"contamination_weight": 2.0' in text and '"consolidator": "metawrap2"' in text


def test_a_second_weight_reuses_every_score_from_the_shared_cache(tmp_path, monkeypatch):
    db = tmp_path / "uniref100.KO.1.dmnd"
    db.write_text("db")
    config = tmp_path / "metawrap2.toml"
    config.write_text(
        '[settings]\nuse_conda_envs = false\nscore_cache_dir = "%s"\n\n[databases]\n'
        'CHECKM2DB = "%s"\n' % (tmp_path / "cache", db)
    )
    monkeypatch.setenv("METAWRAP2_CONFIG", str(config))
    monkeypatch.setattr(scorecache, "tool_version", lambda *a: "checkm2=1.1.0=b")
    calls = []

    def counting(bins_dir, **kwargs):
        calls.append(os.path.basename(bins_dir))
        return _fake_run_checkm(bins_dir, **kwargs)

    monkeypatch.setattr(checkm_mod, "run_checkm", counting)
    monkeypatch.setattr(_common, "_require_env_manager", lambda env: None)

    first = _run(tmp_path / "w5", ["--checkm-jobs", "3"])
    assert calls
    del calls[:]
    second = _run(tmp_path / "w2", ["--contamination-weight", "2"])
    assert calls == [], "the pool is the same, so every candidate was scored by the first run"
    assert _final_bins(first) == _final_bins(second)


def _resolved(extra, capsys=None):
    """What resolve_consolidator makes of *extra*, and the note it printed."""
    args = bin_refinement._parse_args(["-o", "o", "-A", "a", "-B", "b", *extra])
    chosen = bin_refinement.resolve_consolidator(args)
    said = capsys.readouterr().out if capsys is not None else ""
    return chosen, args, " ".join(said.split())


def test_no_flag_gives_the_default_with_no_note(capsys):
    chosen, args, said = _resolved([], capsys)
    assert chosen == "metawrap2" == bin_refinement.DEFAULT_CONSOLIDATOR
    assert args.contamination_weight == CONTAMINATION_WEIGHT
    assert args.refine_combinations == bin_refinement.DEFAULT_COMBINATIONS
    assert "older name" not in said and "applies to" not in said


@pytest.mark.parametrize(
    "extra, chosen, note",
    [
        (["--consolidator", "metawrap2"], "metawrap2", None),
        (["--consolidator", "hybrid"], "metawrap2", "older name; using --consolidator metawrap2"),
        (["--consolidator", "metawrap"], "metawrap1", "older name; using --mw1-refinement"),
        (["--mw1-refinement"], "metawrap1", None),
        (["--consolidator", "metawrap", "--mw1-refinement"], "metawrap1", "older name"),
        (["--consolidator", "dastool"], "dastool", None),
    ],
)
def test_consolidator_names_and_aliases(capsys, extra, chosen, note):
    got, _args, said = _resolved(extra, capsys)
    assert got == chosen
    if note:
        assert note in said
    else:
        assert "older name" not in said


@pytest.mark.parametrize(
    "extra, flag",
    [
        (["--skip-checkm"], "--skip-checkm"),
        (["--skip-refinement"], "--skip-refinement"),
        (["--skip-consolidation"], "--skip-consolidation"),
        (["--refine-combinations", "pairs"], "--refine-combinations"),
        (["--keep-ambiguous"], "--keep-ambiguous"),
        (["--remove-ambiguous"], "--remove-ambiguous"),
    ],
)
def test_an_option_only_mw1_uses_falls_back_to_it_with_a_note(capsys, extra, flag):
    chosen, _args, said = _resolved(extra, capsys)
    assert chosen == bin_refinement.MW1_CONSOLIDATOR
    assert "%s applies to metaWRAP 1's refinement, so this run uses it" % flag in said
    # With --mw1-refinement given, there is nothing to say.
    chosen, _args, said = _resolved(extra + ["--mw1-refinement"], capsys)
    assert chosen == bin_refinement.MW1_CONSOLIDATOR and "applies to" not in said


@pytest.mark.parametrize(
    "extra, message",
    [
        (["--consolidator", "metawrap2", "--skip-checkm"], "metawrap2 and --skip-checkm cannot"),
        (["--consolidator", "hybrid", "--skip-refinement"], "metawrap2 and --skip-refinement"),
        (["--consolidator", "metawrap2", "--keep-ambiguous"], "and --keep-ambiguous cannot"),
        (["--consolidator", "metawrap2", "--refine-combinations", "all"], "--refine-combinations"),
        (["--mw1-refinement", "--consolidator", "dastool"], "--consolidator dastool cannot"),
        (["--mw1-refinement", "--consolidator", "metawrap2"], "--consolidator metawrap2 cannot"),
        (
            ["--mw1-refinement", "--contamination-weight", "2"],
            "--contamination-weight cannot be used with --mw1-refinement",
        ),
        (
            ["--skip-refinement", "--contamination-weight", "2"],
            "--contamination-weight cannot be used with --mw1-refinement",
        ),
        (
            ["--consolidator", "dastool", "--contamination-weight", "2"],
            "--contamination-weight cannot be used with --consolidator dastool",
        ),
    ],
)
def test_conflicting_options_are_refused(capsys, extra, message):
    with pytest.raises(SystemExit):
        _resolved(extra)
    captured = capsys.readouterr()
    assert message in " ".join((captured.out + captured.err).replace("*", " ").split())


def test_the_old_names_are_accepted_but_not_listed(capsys):
    with pytest.raises(SystemExit):
        bin_refinement._parse_args(["-h"])
    text = capsys.readouterr().out
    assert "{metawrap2,dastool}" in text and "--mw1-refinement" in text
    assert "hybrid" not in text.split("--consolidator", 1)[1].split("--mw1-refinement", 1)[0]
    with pytest.raises(SystemExit):
        bin_refinement._parse_args(["-o", "o", "-A", "a", "--consolidator", "nonsense"])
    assert "invalid choice" in capsys.readouterr().err


def test_dastool_still_needs_an_assembly_and_the_default_does_not(tmp_path, scored_dirs, capsys):
    with pytest.raises(SystemExit):
        _run(tmp_path, ["--consolidator", "dastool"])
    captured = capsys.readouterr()
    assert "needs -a/--assembly" in captured.out + captured.err
    args = bin_refinement._parse_args(["-o", "o", "-A", "a"])
    assert args.assembly is None and args.consolidator is None


def test_a_mw1_run_cannot_be_resumed_with_the_default_and_the_reverse(
    tmp_path, scored_dirs, capsys
):
    # A metaWRAP 1 run, including any run from before 2.1.6, is known by its finished
    # Binning_refiner step and writes no consolidator file.
    from argparse import Namespace

    from metawrap2.checkpoint import Checkpoint

    old = str(tmp_path / "old")
    Checkpoint(old, resume=False).done("refine")
    resumed = Checkpoint(old, resume=True)
    with pytest.raises(SystemExit):
        bin_refinement._check_resumed_consolidator(
            old, Namespace(consolidator="metawrap2"), resumed
        )
    captured = capsys.readouterr()
    assert "started with metaWRAP 1's refinement" in captured.out + captured.err
    bin_refinement._check_resumed_consolidator(old, Namespace(consolidator="metawrap1"), resumed)
    assert not os.path.exists(os.path.join(old, "work_files", bin_refinement.CONSOLIDATOR_FILE))

    _run(tmp_path / "new")
    with pytest.raises(SystemExit):
        _run(tmp_path / "new", ["--mw1-refinement"], resume=True)
    captured = capsys.readouterr()
    assert "started with the default consolidator" in captured.out + captured.err
