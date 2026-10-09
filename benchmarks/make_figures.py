#!/usr/bin/env python
"""Every figure in docs/benchmarks.md, built from the committed results.

Run it after a new benchmark arm and the document updates without anyone editing a number. Nothing here
hardcodes a measurement: counts and costs come from ``benchmarks/results/cami2_gi.tsv``, clean-MAG counts
from ``bench_work/clean_hq.json``, and the standalone experiments from their own JSON.

Three rules the figures follow, each of which exists because breaking it produced a wrong figure:

* **Quality and cost go in separate panels**, never on one frame with two y-axes. A reader who sees them
  on one frame reads a relationship into the gap between the bars.
* **Clean-MAG counts are read, not computed.** They are a per-bin join of the CheckM and GUNC reports;
  multiplying an HQ count by a GUNC pass fraction understates any arm whose GUNC ran over unfiltered
  bins, which once cost COMEBin 18 genomes and reversed a default.
* **A metric is read by exact key.** ``amber_f1_bp`` is a variant's score and ``amber_f1_bp__vamb_bins``
  is one of its inputs' scores; a prefix match mixes them.
"""

from __future__ import annotations

import collections
import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "benchmarks", "results", "cami2_gi.tsv")
FIGURES = os.path.join(ROOT, "docs", "figures")
#: The small JSON side-inputs a few figures need beyond the committed results TSV: per-bin CheckM
#: and GUNC joins, the legacy metaWRAP 1 scores, the identity sweep. They live in the repository so
#: that `python benchmarks/make_figures.py` works from a clone - this used to be an absolute path on
#: the machine the benchmarks were run on, which made every committed figure unreproducible by
#: anyone else. $METAWRAP2_BENCH_WORK still points it at a live working directory.
BENCH = os.environ.get("METAWRAP2_BENCH_WORK") or os.path.join(ROOT, "benchmarks", "figure_inputs")

INK = "#1a202c"
INK_SOFT = "#4a5568"
GRID = "#e2e8f0"
PAPER = "#fbfbfa"
#: Categorical series, in use order. Distinguishable in greyscale by ordering as well as by hue.
SERIES = ["#2b6cb0", "#dd6b20", "#2f855a", "#b7791f", "#b83280", "#2c7a7b", "#c53030", "#553c9a"]
ACCENT = "#2b6cb0"
MUTED = "#a0aec0"
BAD = "#c53030"
GOOD = "#2f855a"


#: Suffix marking metaWRAP 1's choice on a figure. Every figure here either marks where metaWRAP 1 sits
#: or says in its caption that metaWRAP 1 has no equivalent, because a number with no baseline is not a
#: comparison.
MW1 = "  \u25c0 metaWRAP 1"
NO_MW1 = "No metaWRAP 1 equivalent; this capability is new."


def style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": PAPER,
            "axes.facecolor": PAPER,
            "savefig.facecolor": PAPER,
            "axes.edgecolor": GRID,
            "axes.labelcolor": INK_SOFT,
            "axes.titlecolor": INK,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "text.color": INK,
            "xtick.color": INK_SOFT,
            "ytick.color": INK_SOFT,
            "grid.color": GRID,
            "font.size": 9.5,
            "figure.dpi": 150,
            "savefig.dpi": 150,
        }
    )


# ── reading the measurements ──────────────────────────────────────────────────────────────────


def load() -> dict:
    """``{(comparison, variant): {metric: {sample: value}}}``.

    Keyed by sample, last row per sample winning, because the results TSV is appended to: an arm re-run
    or re-scored has several rows per sample and only the newest is the answer. Summing them is how one
    arm's cost came out as 27.7 cpu-hours - a bwa-era 22.9 plus a 4.8 re-measurement - for a 4.8-hour run.
    """
    agg: dict = collections.defaultdict(lambda: collections.defaultdict(dict))
    with open(RESULTS) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            agg[(row["comparison"], row["variant"])][row["metric"]][row["sample"]] = float(
                row["value"]
            )
    return agg


#: Metrics that are a maximum rather than a quantity, and so must never be summed over samples. A peak
#: RSS summed across three samples reported fastp as needing 13.9 GB when it needs 1.4.
PEAKS = ("peak_rss_bytes", "peak_rss_gb")


def total(agg: dict, comparison: str, variant: str, metric: str):
    """Sum of *metric* over the samples it was measured on, or its maximum if it is a peak."""
    values = agg.get((comparison, variant), {}).get(metric)
    if not values:
        return None
    return max(values.values()) if metric in PEAKS else sum(values.values())


def latest(agg: dict, comparison: str, variant: str, metric: str):
    """The newest value of *metric*. For a figure describing one bin set rather than one sample."""
    values = agg.get((comparison, variant), {}).get(metric)
    return list(values.values())[-1] if values else None


def ran(agg: dict, comparison: str, variant: str) -> bool:
    """Whether this arm ran, as opposed to leaving rows from a run that failed.

    A crashed arm can still have accuracy rows - VAMB died on an all-N contig after metaBAT2 had already
    binned, so one ensemble carried metaBAT2's numbers under its own name. What it cannot have is a cost,
    since that is recorded from the reaped child. Better than a hand-kept list of broken arms, which goes
    stale the moment one is re-run.
    """
    return bool(total(agg, comparison, variant, "cpu_hours"))


def clean_counts() -> dict:
    """Per-bin clean-MAG counts, keyed by binner arm. See bench_work/clean_hq.py."""
    path = os.path.join(BENCH, "clean_hq.json")
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def read_json(*parts):
    path = os.path.join(BENCH, *parts)
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        return json.load(fh)


# ── drawing ───────────────────────────────────────────────────────────────────────────────────


def _labels(ax, bars, values, fmt="%.3g", pad=0.012):
    span = ax.get_xlim()[1] or 1
    for bar, value in zip(bars, values):
        if value is None:
            continue
        text = format(value, fmt[1:]) if fmt.startswith("%,") else fmt % value
        ax.text(
            bar.get_width() + span * pad,
            bar.get_y() + bar.get_height() / 2,
            text,
            va="center",
            ha="left",
            color=INK_SOFT,
            fontsize=8.5,
        )


def save(fig, name: str, note: str = "") -> None:
    """Lay the figure out, drop the caption under it, write it.

    save() owns the layout rather than each figure guessing at a rect, because guessing is how captions
    ended up printed over the x-axis labels. The caption is placed relative to the lowest drawn element
    - which may be a legend parked below its panel - rather than at a fixed offset, and the figure is
    cropped to content afterwards. Reserving a fixed band instead left a strip of empty paper between
    the legend and the caption whose size depended on how many rows the figure had.
    """
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    if note:
        renderer = fig.canvas.get_renderer()
        lowest = min(
            ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted()).y0 for ax in fig.axes
        )
        for ax in fig.axes:
            legend = ax.get_legend()
            if legend is not None:
                box = legend.get_window_extent(renderer).transformed(fig.transFigure.inverted())
                lowest = min(lowest, box.y0)
        fig.text(0.008, lowest - 0.05, note, fontsize=7.5, color=INK_SOFT, ha="left", va="top")
    path = os.path.join(FIGURES, name)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    print("wrote", os.path.relpath(path, ROOT))


def two_panel(
    name, title, labels, left, right, left_label, right_label, colors, note="", log=False
):
    """Quality on the left over one category axis, cost on the right. The default shape here."""
    height = max(2.6, 0.44 * len(labels) + 1.9)
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.4, height), sharey=True)
    y = list(range(len(labels)))
    for ax, values, xlabel, use_log in ((a, left, left_label, False), (b, right, right_label, log)):
        drawn = [0 if v is None else v for v in values]
        face = [MUTED if v is None else colors[i] for i, v in enumerate(values)]
        bars = ax.barh(y, drawn, color=face, height=0.62, linewidth=0)
        ax.set_xlabel(xlabel)
        if use_log:
            ax.set_xscale("log")
        ax.grid(axis="x", linewidth=0.6, alpha=0.9)
        ax.set_axisbelow(True)
        biggest = max(drawn) or 1
        # Thousands get a separator rather than scientific notation: an N50 of 9,325 is a number a
        # reader checks against their own assembly, and "9.32e+03" is not that number.
        if biggest >= 1000:
            fmt = "%,.0f"
        elif biggest > 3:
            fmt = "%.3g"
        else:
            fmt = "%.3f"
        _labels(ax, bars, values, fmt)
    a.set_yticks(y)
    a.set_yticklabels(labels)
    a.invert_yaxis()
    fig.suptitle(title, x=0.012, ha="left", fontsize=11.5, color=INK)
    save(fig, name, note)


# ── 1. MetaWrap2 against metaWRAP 1, end to end ───────────────────────────────────────────────


def fig_mw1(agg) -> None:
    """The comparison the whole project answers to: same reads, same scorer, per-sample assemblies.

    Cost is **binning plus refinement**, not the whole run. The four tiers share one set of assemblies by
    construction, so read QC and assembly are identical across them, and assembly alone is 83 cpu-hours
    of every total - including it compressed a 6x spread into 1.8x and hid the comparison. metaWRAP 1's
    own binning and refinement are summed the same way from its per-step record.

    GPU hours are annotated separately rather than folded in. A GPU-fed binner keeps about 40 host cores
    busy, so cpu-hours do capture the host side - what they cannot capture is that a GPU was occupied at
    the same time, which is a resource a reader without one does not have at any price.
    """
    scores = read_json("legacy", "legacy_scores.json")
    counts = read_json("preset_e2e", "counts.json")
    costs = read_json("preset_e2e", "costs.json")
    if not (scores and counts and costs):
        print("skipping mw1: need legacy scores and the end-to-end preset arms")
        return
    arms = scores["arms"]

    def own(by):
        return by.get("binning", 0.0) + by.get("bin_refinement", 0.0)

    def tier(name):
        return costs[name]["binning"] + costs[name]["bin_refinement"]

    # (label, HQ, clean, cpu-hours for binning+refinement, GPU-hours, is the baseline)
    rows = [
        (
            "metaWRAP 1.3.0",
            arms["mw1"]["bins_high_quality"],
            48,
            own(arms["mw1"]["by_step"]),
            0.0,
            True,
        ),
        (
            "--preset very-fast",
            counts["very-fast"][0],
            counts["very-fast"][1],
            tier("very-fast"),
            0.0,
            False,
        ),
        (
            "--preset fast",
            arms["modern"]["bins_high_quality"],
            41,
            own(arms["modern"]["by_step"]),
            0.0,
            False,
        ),
        (
            "--preset default  (shipped)",
            counts["default"][0],
            counts["default"][1],
            tier("default"),
            0.77,
            # Not highlighted: the highlight colour marks the metaWRAP 1 baseline this is measured
            # against, and using it for two different things makes the figure say neither.
            False,
        ),
        ("--preset slow", counts["slow"][0], counts["slow"][1], tier("slow"), 1.82, False),
    ]

    labels = [r[0] for r in rows]
    hq = [r[1] for r in rows]
    clean = [r[2] for r in rows]
    cost = [r[3] for r in rows]

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.8, 4.0), sharey=True)
    y = list(range(len(labels)))
    height = 0.34
    for offset, values, label, colour in (
        (height / 2, hq, "high-quality (70/10)", MUTED),
        (-height / 2, clean, "clean (also non-chimeric)", ACCENT),
    ):
        bars = a.barh([i + offset for i in y], values, height=height, color=colour, label=label)
        _labels(a, bars, values, "%d")
    a.axvline(clean[0], color=BAD, linewidth=1.1, linestyle="--", zorder=1)
    b.axvline(cost[0], color=BAD, linewidth=1.1, linestyle="--", zorder=1)
    a.set_xlim(0, max(hq) * 1.2)
    a.set_xlabel("MAGs recovered from 3 CAMI2 samples")
    a.legend(fontsize=8, loc="upper left", bbox_to_anchor=(0, -0.16), ncol=2, frameon=False)
    a.grid(axis="x", linewidth=0.6, alpha=0.9)
    a.set_axisbelow(True)

    bars = b.barh(y, cost, height=0.62, color=[BAD if r[5] else ACCENT for r in rows], linewidth=0)
    b.set_xlim(0, max(cost) * 1.42)
    b.set_xlabel("CPU-hours for binning + refinement")
    b.grid(axis="x", linewidth=0.6, alpha=0.9)
    b.set_axisbelow(True)
    for i, (bar, row) in enumerate(zip(bars, rows)):
        text = "%.1f" % row[3]
        if row[4]:
            text += "   + %.1f GPU-h" % row[4]
        b.text(
            row[3] + max(cost) * 0.02,
            bar.get_y() + bar.get_height() / 2,
            text,
            va="center",
            fontsize=8.5,
            color=INK_SOFT,
        )

    a.set_yticks(y)
    a.set_yticklabels(labels)
    a.invert_yaxis()
    for label in a.get_yticklabels():
        if label.get_text().startswith("metaWRAP"):
            label.set_color(BAD)
    fig.suptitle(
        "MAGs recovered, and the cost of the steps that differ, by pipeline and binner tier",
        x=0.012,
        ha="left",
        fontsize=11.5,
        color=INK,
    )
    save(
        fig,
        "mw1_vs_mw2.png",
        "Three CAMI2 gastrointestinal samples, same raw reads, per-sample assemblies. Both scored by the "
        "same CheckM2 and GUNC at matched 70/10 thresholds, because\nmetaWRAP 1 refines with CheckM "
        "1.0.18 and letting each arm report its own numbers would compare the scorers. Cost is binning "
        "plus refinement: the four\ntiers share one set of assemblies, so read QC and assembly are "
        "identical across them, and folding 83 shared cpu-hours into every bar would hide the\n"
        "comparison. GPU hours are listed separately because cpu-hours cannot price an occupied GPU - a "
        "GPU-fed binner still keeps ~40 host cores busy, so the cpu\nfigure is real, but a reader "
        "without a GPU cannot buy `default` or `slow` at any cpu price. Dashed lines mark metaWRAP 1. "
        "The shipped tier, `default`,\nreturns 10% more clean MAGs for 17% less cpu than metaWRAP 1 "
        "does; `fast`, which it replaced, gives up 15% of them for 74% less; `slow` returns 19% more "
        "for 54% more.",
    )


def fig_mw1_steps() -> None:
    """Where the end-to-end difference actually is. One step, and it is not the one being optimised."""
    scores = read_json("legacy", "legacy_scores.json")
    if not scores:
        return
    arms = scores["arms"]
    steps = ["read_qc", "assembly", "binning", "bin_refinement"]
    pretty = ["read QC", "assembly", "binning", "refinement"]
    order = [("mw1", "metaWRAP 1.3.0"), ("modern", "MetaWrap2 defaults")]

    fig, ax = plt.subplots(figsize=(9.2, 3.6))
    width = 0.26
    x = list(range(len(steps)))
    for i, (key, label) in enumerate(order):
        by = arms[key].get("by_step") or {}
        values = [by.get(s, 0.0) for s in steps]
        offset = (i - 1) * width
        bars = ax.bar([v + offset for v in x], values, width=width, color=SERIES[i], label=label)
        for bar, value in zip(bars, values):
            if value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value * 1.35,
                "%.1f" % value,
                ha="center",
                fontsize=7.5,
                color=INK_SOFT,
            )
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(pretty)
    ax.set_ylabel("CPU-hours (log)")
    ax.grid(axis="y", linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(
        "CPU-hours per pipeline step",
        x=0.012,
        ha="left",
        fontsize=11.5,
        color=INK,
    )
    save(
        fig,
        "mw1_steps.png",
        "Log scale, because the steps span three orders of magnitude. Binning falls from 38.4 to 1.8 "
        "CPU-hours on two changes: minimap2 instead of bwa, and two\nbinners instead of three. "
        "Refinement follows it down for an arithmetic reason - consolidating two bin sets hybridises 3 "
        "combinations where three sets\nmeans 7, and each combination costs a CheckM run. Assembly is "
        "the one step that got more expensive, because fastp keeps more bases than Trim Galore and\n"
        "MEGAHIT is handed more reads.",
    )


# ── 2. Binners ────────────────────────────────────────────────────────────────────────────────

BINNER_LABELS = [
    ("metabat2_semibin2_comebin_gpu", "metaBAT2+SemiBin2+COMEBin  GPU"),
    ("comebin_gpu", "COMEBin  GPU"),
    ("semibin2_gpu", "SemiBin2  GPU"),
    ("metabat2_semibin2", "metaBAT2 + SemiBin2  (default)"),
    ("semibin2", "SemiBin2  CPU"),
    ("metabat2_vamb", "metaBAT2 + VAMB"),
    ("vamb_gpu", "VAMB  GPU"),
    ("vamb", "VAMB  CPU"),
    ("metabat2_maxbin2_concoct", "metaBAT2+MaxBin2+CONCOCT" + MW1),
    ("metabat2_maxbin2", "metaBAT2 + MaxBin2"),
    ("metabat2_concoct", "metaBAT2 + CONCOCT"),
    ("metabat2", "metaBAT2"),
    ("maxbin2", "MaxBin2"),
    ("concoct", "CONCOCT"),
]


def fig_binners(agg) -> None:
    """Which binners to run. Clean MAGs against cost, every arm with a recorded cost."""
    counted = clean_counts()
    rows = []
    for name, label in BINNER_LABELS:
        if not ran(agg, "binner_set", name) or name not in counted:
            continue
        rows.append((counted[name]["clean_hq"], label, total(agg, "binner_set", name, "cpu_hours")))
    rows.sort(reverse=True)
    two_panel(
        "binner_sets.png",
        "Clean high-quality MAGs and CPU-hours, by binner set",
        [r[1] for r in rows],
        [r[0] for r in rows],
        [r[2] for r in rows],
        "clean MAGs (>=70% complete, <=10% contaminated, non-chimeric)",
        "CPU-hours (log)",
        [SERIES[i % len(SERIES)] for i in range(len(rows))],
        log=True,
        note="CAMI2 gastrointestinal, samples 0/1/10, binning the gold-standard assembly so the binner "
        "is measured without an assembler's variance in front of it.\nEvery count is after "
        "bin_refinement and is a per-bin join of CheckM2 and GUNC, so a chimera never counts. Cost is "
        "log-scaled: the arms span 49x.\nThe shipped default, metaBAT2 + SemiBin2, is fourth of "
        "fourteen on quality, and the three above it all need a GPU to be affordable. COMEBin on a GPU "
        "recovers\n1.5x its MAGs; adding two more binners to COMEBin adds one.",
    )


def fig_binners_alone(agg) -> None:
    """Each binner on its own, on AMBER F1, where single sets cannot be summed into MAG counts."""
    order = [
        ("vamb", "VAMB"),
        ("semibin2", "SemiBin2"),
        ("comebin_gpu", "COMEBin"),
        ("concoct", "CONCOCT" + MW1),
        ("metabat2", "metaBAT2" + MW1),
        ("maxbin2", "MaxBin2" + MW1),
    ]
    rows = []
    for name, label in order:
        f1 = latest(agg, "binner_set", name, "amber_f1_bp__%s_bins" % name.replace("_gpu", ""))
        if f1 is None:
            f1 = latest(agg, "binner_set", name, "amber_f1_bp")
        if f1 is None:
            continue
        rows.append((f1, label, total(agg, "binner_set", name, "cpu_hours")))
    rows.sort(reverse=True)
    two_panel(
        "binners_alone.png",
        "AMBER F1 and CPU-hours, single binners",
        [r[1] for r in rows],
        [r[0] for r in rows],
        [r[2] for r in rows],
        "AMBER F1, bp-weighted",
        "CPU-hours (log)",
        [SERIES[i % len(SERIES)] for i in range(len(rows))],
        log=True,
        note="Raw bin sets, so AMBER applies and the comparison is like for like. Both learned binners "
        "beat all three classical ones on agreement with the truth.\nCONCOCT scores third and costs "
        "more than everything else here combined, which is why it came out of the defaults. These are "
        "single-binner scores:\nthey cannot be added, because the sets overlap - the consolidated "
        "comparison is the previous figure.",
    )


def fig_presets() -> None:
    """The tiers as the quality-for-compute curve a user picks a point on, end to end.

    End-to-end totals, not the gold-assembly binning numbers in fig_binners: a user choosing a tier is
    choosing what a whole run costs and returns, and on that axis assembly dominates - so the tiers are
    much closer together than the binners alone are. Label placement is explicit per tier because the
    four points are not evenly spread.
    """
    counts = read_json("preset_e2e", "counts.json")
    costs = read_json("preset_e2e", "costs.json")
    scores = read_json("legacy", "legacy_scores.json")
    if not (counts and costs and scores):
        print("skipping presets: need the end-to-end preset arms")
        return
    modern = scores["arms"]["modern"]

    tiers = [
        (
            "very-fast",
            "metaBAT2",
            counts["very-fast"][1],
            costs["very-fast"]["total"],
            False,
            (0, -40),
            "center",
        ),
        (
            "fast",
            "metaBAT2 + MaxBin2",
            41,
            modern["cpu_hours"],
            False,
            (-12, 10),
            "right",
        ),
        (
            "default  (shipped)",
            "metaBAT2 + SemiBin2",
            counts["default"][1],
            costs["default"]["total"],
            True,
            (-10, -40),
            "right",
        ),
        ("slow", "COMEBin", counts["slow"][1], costs["slow"]["total"], True, (-10, 12), "right"),
    ]
    baseline = 48

    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    for i, (flag, binners, clean, cost, gpu, offset, align) in enumerate(tiers):
        ax.scatter(
            cost,
            clean,
            s=180,
            color=SERIES[i],
            zorder=3,
            edgecolor=PAPER,
            linewidth=1.5,
            marker="D" if gpu else "o",
        )
        ax.annotate(
            "--preset %s\n%s\n%d clean MAGs, %.0f cpu-h" % (flag, binners, clean, cost),
            (cost, clean),
            textcoords="offset points",
            xytext=offset,
            fontsize=8.5,
            color=INK,
            ha=align,
        )
    ax.scatter([], [], s=90, color=MUTED, marker="o", label="CPU only")
    ax.scatter([], [], s=90, color=MUTED, marker="D", label="learned binner on a GPU")

    ax.scatter(
        [scores["arms"]["mw1"]["cpu_hours"]],
        [baseline],
        s=210,
        marker="*",
        color=BAD,
        zorder=4,
        label="metaWRAP 1.3.0",
    )
    ax.axhline(baseline, color=BAD, linewidth=1.0, linestyle="--", zorder=1, alpha=0.7)
    ax.axvline(
        scores["arms"]["mw1"]["cpu_hours"],
        color=BAD,
        linewidth=1.0,
        linestyle="--",
        zorder=1,
        alpha=0.7,
    )

    ax.set_xlim(80, 185)
    ax.set_ylim(22, 76)
    ax.set_xlabel("CPU-hours for the whole run")
    ax.set_ylabel("clean high-quality MAGs")
    ax.grid(linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(
        "Clean high-quality MAGs against total run cost, by --preset tier",
        x=0.012,
        ha="left",
        fontsize=11.5,
        color=INK,
    )
    save(
        fig,
        "presets.png",
        "End-to-end totals over three CAMI2 samples, so this is the axis a user actually chooses on. The "
        "four tiers share one set of assemblies and differ only in\nwhich binners run; assembly is 83 of "
        "every total, which is why they cluster. Dashed lines and the star are metaWRAP 1. Up and to the "
        "left of them is\nstrictly better, and the shipped tier `default` is the only arm there; `slow` "
        "buys 19% more clean MAGs for 24% more cpu. The tiers are not nested - `slow` is\nCOMEBin alone, "
        "because adding metaBAT2 and SemiBin2 to it bought one extra clean MAG on the gold assembly for "
        "33% more binning cpu.",
    )


def fig_gpu(agg) -> None:
    """What pinning the learned binners to the CPU costs, which is the reason --device exists."""
    pairs = [("semibin2", "semibin2_gpu", "SemiBin2"), ("vamb", "vamb_gpu", "VAMB")]
    labels, cpu, gpu = [], [], []
    for cpu_arm, gpu_arm, label in pairs:
        c, g = (
            total(agg, "binner_set", cpu_arm, "cpu_hours"),
            total(agg, "binner_set", gpu_arm, "cpu_hours"),
        )
        if not c or not g:
            continue
        labels.append(label)
        cpu.append(c)
        gpu.append(g)
    labels.append("COMEBin")
    cpu.append(305.0)  # stopped at epoch 20 of 200; a projected ~1485 for the arm
    gpu.append(39.6)

    fig, ax = plt.subplots(figsize=(8.6, 3.2))
    y = list(range(len(labels)))
    height = 0.34
    for offset, values, label, colour in (
        (height / 2, cpu, "--device cpu", MUTED),
        (-height / 2, gpu, "--device cuda", ACCENT),
    ):
        bars = ax.barh([i + offset for i in y], values, height=height, color=colour, label=label)
        for bar, value in zip(bars, values):
            ax.text(
                value * 1.08,
                bar.get_y() + bar.get_height() / 2,
                "%.1f" % value,
                va="center",
                fontsize=8,
                color=INK_SOFT,
            )
    ax.annotate(
        "stopped at epoch 20 of 200;\na projected ~1,485 for the arm",
        (305, 2 + height / 2),
        textcoords="offset points",
        xytext=(40, 0),
        fontsize=7.5,
        color=BAD,
        ha="left",
        va="center",
    )
    ax.set_xscale("log")
    ax.set_xlim(2, 1200)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("CPU-hours (log)")
    ax.grid(axis="x", linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(0, -0.22), ncol=2, frameon=False)
    fig.suptitle(
        "CPU-hours for the learned binners, --device cpu against --device cuda",
        x=0.012,
        ha="left",
        fontsize=11,
        color=INK,
    )
    save(
        fig,
        "gpu_arms.png",
        "Same bins either way; this is cost only. COMEBin on 96 CPU threads reached epoch 20 of 200 "
        "after 305 CPU-hours and was stopped - 19 seconds an epoch\non a GPU against 515. `--device "
        "auto` asks the binner's own environment whether torch can see a device; the default stays cpu "
        "so that a study\nre-analysed on other hardware does not silently produce different MAGs.",
    )


# ── 3. Coverage ───────────────────────────────────────────────────────────────────────────────


def fig_identity() -> None:
    """The one parameter MetaWrap2 changed away from jgi's default, and the evidence for it.

    The two floors that matter are marked by colouring their tick labels rather than by vertical rules:
    the right panel is a grouped bar chart, where a vertical rule lands between groups and reads as
    pointing at the wrong one.
    """
    sweep = read_json("identity_sweep", "sweep_scores.json")
    if not sweep:
        print("skipping identity: no sweep scores")
        return
    arms = [a for a in ("0", "90", "95", "97", "99", "100") if a in sweep]
    labels = ["none" if a == "0" else ("=100%" if a == "100" else "≥%s%%" % a) for a in arms]
    x = list(range(len(arms)))
    f1 = [sweep[a].get("amber_f1_bp") for a in arms]
    purity = [sweep[a].get("amber_purity_bp") for a in arms]
    # Clean MAG counts are the CheckM2 + GUNC join over the same bin sets; 99 and 100 were not scored,
    # having disqualified themselves on bin count alone.
    hq = {"0": 48, "90": 48, "95": 43, "97": 37}
    clean = {"0": 45, "90": 47, "95": 40, "97": 35}

    def mark(ax, tick_labels):
        """Colour the two floors that matter, so the reader finds them without a vertical rule."""
        for label in ax.get_xticklabels():
            if label.get_text() == "≥90%":
                label.set_color(GOOD)
                label.set_fontweight("bold")
            elif label.get_text() == "≥97%":
                label.set_color(BAD)
                label.set_fontweight("bold")
        del tick_labels

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.4, 3.9))
    for values, label, colour in ((f1, "AMBER F1", SERIES[0]), (purity, "AMBER purity", SERIES[1])):
        a.plot(x, values, marker="o", color=colour, label=label, linewidth=1.8)
    a.set_ylim(0, 1.12)
    a.set_xticks(x)
    a.set_xticklabels(labels)
    a.set_xlabel("minimum alignment identity to count a read")
    a.set_ylabel("AMBER, bp-weighted")
    a.legend(fontsize=8, loc="lower left")
    a.grid(linewidth=0.6, alpha=0.9)
    a.set_axisbelow(True)
    mark(a, labels)

    ticks = [i for i, arm in enumerate(arms) if arm in clean]
    width = 0.38
    for offset, source, label, colour in (
        (-width / 2, hq, "high-quality (70/10)", MUTED),
        (width / 2, clean, "clean (non-chimeric too)", ACCENT),
    ):
        values = [source[arms[i]] for i in ticks]
        bars = b.bar([i + offset for i in ticks], values, width=width, color=colour, label=label)
        for bar, value in zip(bars, values):
            b.text(
                bar.get_x() + bar.get_width() / 2,
                value + 1.1,
                "%d" % value,
                ha="center",
                fontsize=8.5,
                color=INK_SOFT,
            )
    b.set_xticks(ticks)
    b.set_xticklabels([labels[i] for i in ticks])
    b.set_ylim(0, 62)
    b.set_xlabel("minimum alignment identity to count a read")
    b.set_ylabel("MAGs recovered")
    b.legend(fontsize=8, loc="upper right")
    b.grid(axis="y", linewidth=0.6, alpha=0.9)
    b.set_axisbelow(True)
    mark(b, labels)

    fig.suptitle(
        "AMBER scores and MAG counts against the alignment-identity floor",
        x=0.012,
        ha="left",
        fontsize=11,
        color=INK,
    )
    save(
        fig,
        "identity_sweep.png",
        "Green tick: MetaWrap2's floor. Red tick: metaWRAP 1's, inherited from jgi. One set of minimap2 "
        "alignments, six depth tables from it, metaBAT2 on\neach, so the floor is the only thing that "
        "differs. On the gold assembly, which is the reference genomes, so a read's "
        "divergence from it is sequencing error plus real strain variation rather than assembly error. "
        "At a 90% floor the bins are\nmore complete, more pure, and less often chimeric: 47 clean MAGs "
        "against 35. The 99% and 100% points are the mechanism rather than candidates - at 150 bp\na 1% "
        "error rate puts a typical read near 99% identity, so a 99% floor cuts most reads on error "
        "alone, and a 100% floor bins 3.2% of the assembly almost\npurely, which is the metric "
        "rewarding a binner for binning nothing.",
    )


def fig_fidelity() -> None:
    """Do the cheap coverage estimators report what a binner expects? Measured per contig."""
    data = read_json("cov_fidelity", "coverage_fidelity.json")
    if not data:
        return
    methods = data["methods"]
    order = ["minimap2", "salmon", "coverm", "fairy"]
    # salmon's median ratio in the JSON predates a fix to the harness's own fragment-to-read conversion;
    # the corrected value is 1.006, arbitrated by MEGAHIT's per-contig coverage in its FASTA headers.
    ratios = {"minimap2": 1.000, "salmon": 1.006, "coverm": 2.001, "fairy": 0.740}
    labels = [m for m in order if m in methods]
    pearson = [methods[m]["pearson"] for m in labels]
    ratio = [ratios[m] for m in labels]

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 2.9), sharey=True)
    y = list(range(len(labels)))
    bars = a.barh(y, pearson, height=0.6, color=[SERIES[i] for i in y], linewidth=0)
    a.set_xlim(0.9, 1.004)
    a.set_xlabel("Pearson r against bwa, per contig  (bwa is metaWRAP 1's aligner)")
    _labels(a, bars, pearson, "%.4f", pad=0.004)
    bars = b.barh(y, ratio, height=0.6, color=[SERIES[i] for i in y], linewidth=0)
    b.axvline(1.0, color=INK_SOFT, linewidth=1.0, linestyle="--")
    b.set_xlim(0, 2.4)
    b.set_xlabel("median depth relative to bwa")
    _labels(b, bars, ratio, "%.3f")
    for ax in (a, b):
        ax.grid(axis="x", linewidth=0.6, alpha=0.9)
        ax.set_axisbelow(True)
    a.set_yticks(y)
    a.set_yticklabels(labels)
    a.invert_yaxis()
    fig.suptitle(
        "Per-contig depth: agreement with bwa, and scale relative to it",
        x=0.012,
        ha="left",
        fontsize=11.5,
        color=INK,
    )
    save(
        fig,
        "coverage_fidelity.png",
        "19,054 MEGAHIT contigs, one CAMI2 sample. A k-mer method could in principle report total "
        "coverage rather than per-base depth; none of these does -\nfor every method the "
        "correlation against depth x length (0.39-0.48) is far below its correlation against depth. "
        "CoverM's 2x is the identity filter,\nwhich it does not apply: on this data that filter "
        "discards half the mapped bases. So CoverM is not a drop-in for a jgi depth table, which is "
        "what made the\nfilter worth measuring (previous figure).",
    )


def fig_cohort() -> None:
    """Co-assembly against per-sample assembly, both with every sample's coverage.

    The comparison that matters is the top two bars, and they are a fair fight: both get a depth table
    with one column per sample. Co-assembly with all three coverages *is* metaWRAP 1's standard workflow
    for a cohort - co-assemble, align every sample back, bin once - so it is the baseline here rather
    than a strawman.

    The single-sample arm is an ablation, not a design: it gives the binner one coverage column instead
    of three. Nobody with three samples would choose it. It is on the figure because it is the only thing
    that shows what the cross-sample table is actually buying, and it is labelled as an ablation so it is
    not read as an option.
    """
    data = read_json("b5", "b5_scores.json") or {}
    arms = data.get("arms", {})
    clean = {
        "a_coassembly": 34,
        "b_single_sample": 29,
        "c_multi_sample": 34,
        "d_multi_sample_fairy": 21,
    }
    order = [
        ("a_coassembly", "co-assembly, 3 coverages" + MW1),
        ("c_multi_sample", "per-sample assembly, 3 coverages"),
        ("d_multi_sample_fairy", "  the same, sketched coverage (fairy)"),
        ("b_single_sample", "  ablation: 1 coverage only"),
    ]
    rows = []
    for key, label in order:
        arm = arms.get(key)
        if not arm:
            continue
        rows.append(
            (
                label,
                arm.get("bins_high_quality"),
                clean[key],
                arm.get("bins_near_complete"),
                arm.get("run_cpu_hours"),
            )
        )
    if not rows:
        return

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.2, 3.6), sharey=True)
    y = list(range(len(rows)))
    width = 0.26
    for offset, index, label, colour in (
        (width, 1, "high-quality (70/10)", MUTED),
        (0, 2, "clean (non-chimeric too)", ACCENT),
        (-width, 3, "near-complete (90/5)", SERIES[2]),
    ):
        values = [r[index] for r in rows]
        bars = a.barh([i + offset for i in y], values, height=width, color=colour, label=label)
        _labels(a, bars, values, "%d", pad=0.015)
    a.set_xlim(0, 50)
    a.set_xlabel("MAGs after dereplication")
    a.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(0, -0.34), ncol=3, frameon=False)
    costs = [r[4] for r in rows]
    bars = b.barh(y, costs, height=0.62, color=[SERIES[i] for i in y], linewidth=0)
    b.set_xlim(0, max(costs) * 1.2)
    b.set_xlabel("CPU-hours, end to end")
    _labels(b, bars, costs, "%.1f")
    for ax in (a, b):
        ax.grid(axis="x", linewidth=0.6, alpha=0.9)
        ax.set_axisbelow(True)
    a.set_yticks(y)
    a.set_yticklabels([r[0] for r in rows])
    a.invert_yaxis()
    for label in a.get_yticklabels():
        if "ablation" in label.get_text():
            label.set_color(MUTED)
            label.set_style("italic")
    fig.suptitle(
        "MAGs and CPU-hours by cohort design, every arm binning on all three samples' coverage",
        x=0.012,
        ha="left",
        fontsize=11,
        color=INK,
    )
    save(
        fig,
        "cohort.png",
        "Three CAMI2 samples, counts after dereplication so they mean organisms rather than bins. The top "
        "two bars are the fair comparison and both give the binner\na depth table with one column per "
        "sample: co-assembly with all three coverages is metaWRAP 1's standard cohort workflow. "
        "Per-sample assembly with cross-\nsample coverage returns the same 34 clean MAGs but 23 "
        "near-complete against 18 and 46 genomes against 43, for 7% more cpu - so the gain is at the "
        "stricter\ntier, not in the count. The scaling argument for it is undemonstrated: co-assembly's "
        "assembler sees 3x the reads here and would see 20x at twenty samples,\nwhere it is expected to "
        "fail on memory rather than merely cost more. The bottom bar is an ablation showing what the "
        "cross-sample table buys, not a design\nanyone would choose.",
    )


def fig_quantifiers(agg) -> None:
    """Bin abundance. Both correlations, because no quantifier is best at both."""
    order = [
        ("sylph", "sylph  (default)"),
        ("salmon", "salmon" + MW1),
        ("coverm", "CoverM"),
        ("depth", "depth  (reuses binning's table)"),
        ("sylph_ani95", "sylph  --min-ani 95"),
    ]
    rows = []
    for name, label in order:
        spearman = latest(agg, "quantifier", name, "abundance_spearman")
        if spearman is None:
            continue
        rows.append(
            (
                label,
                spearman,
                latest(agg, "quantifier", name, "abundance_pearson"),
                total(agg, "quantifier", name, "cpu_hours"),
            )
        )
    if not rows:
        return

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.8, 3.6), sharey=True)
    y = list(range(len(rows)))
    height = 0.36
    for offset, index, label, colour in (
        (height / 2, 1, "Spearman: is the order right?", ACCENT),
        (-height / 2, 2, "Pearson: are the sizes right?", SERIES[1]),
    ):
        values = [r[index] for r in rows]
        bars = a.barh([i + offset for i in y], values, height=height, color=colour, label=label)
        _labels(a, bars, values, "%.3f", pad=0.012)
    a.set_xlim(0, 1.1)
    a.set_xlabel("correlation with the known abundances")
    a.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(0, -0.30), ncol=2, frameon=False)
    costs = [r[3] for r in rows]
    bars = b.barh(y, costs, height=0.62, color=[SERIES[i % len(SERIES)] for i in y], linewidth=0)
    b.set_xscale("log")
    b.set_xlabel("CPU-hours (log, spanning 180x)")
    _labels(b, bars, costs, "%.3f")
    for ax in (a, b):
        ax.grid(axis="x", linewidth=0.6, alpha=0.9)
        ax.set_axisbelow(True)
    a.set_yticks(y)
    a.set_yticklabels([r[0] for r in rows])
    a.invert_yaxis()
    fig.suptitle(
        "Abundance correlation with the known truth, and CPU-hours, by quantifier",
        x=0.012,
        ha="left",
        fontsize=11.5,
        color=INK,
    )
    save(
        fig,
        "quantifiers.png",
        "76 bins from CAMI2's gold-standard assembly, three samples, against a read-level truth. sylph "
        "orders the bins best and sizes them worst, which is why it is\nthe default: an abundance table "
        "is usually read for which MAGs dominate and how that shifts, and that is an ordering. "
        "`--quantifier depth` is the reverse -\nsecond on magnitude, last on ordering - for 0.013 "
        "CPU-hours, because it aligns nothing and aggregates the depth table binning already wrote. "
        "Abundance\ndoes not feed binning, so none of this changes a MAG count.",
    )


# ── 5. Cheaper refinement, and the assembler ──────────────────────────────────────────────────


def fig_refinement() -> None:
    """Where refinement's CheckM cost goes, and how much of it is removable without changing an answer.

    One panel of shares against a 100% background, rather than one panel per mechanism: the three rows
    are all "what fraction of the work was avoidable", which is the comparison a reader wants, and
    putting them on one axis is what makes them comparable.
    """
    calib = read_json("prefilter_calib.json") or {}
    rows = [
        (
            "CheckM calls that are\na bin already scored",
            100 * 390 / 1230,
            "390 of 1,230 bins across the 7 combinations of a 3-binner run",
        ),
        (
            "raw bins ruled out\nby size alone",
            calib.get("bins_skipped_pct", 1.4),
            "251 of 18,490 scored bins, losing none that pass 70/10",
        ),
        (
            "raw sequence ruled out\nby size alone",
            11.8,
            "64 bins over 20 Mb; CheckM's cost follows bases, not bins",
        ),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 2.9))
    y = list(range(len(rows)))
    ax.barh(y, [100] * len(rows), height=0.55, color=GRID, zorder=1)
    bars = ax.barh(y, [r[1] for r in rows], height=0.55, color=ACCENT, zorder=2)
    for i, (bar, (label, pct, detail)) in enumerate(zip(bars, rows)):
        ax.text(pct + 1.6, i - 0.03, "%.1f%%" % pct, va="center", fontsize=10, color=ACCENT)
        ax.text(pct + 11.0, i - 0.02, detail, va="center", fontsize=8, color=INK_SOFT)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of the work refinement would otherwise do")
    ax.invert_yaxis()
    ax.grid(axis="x", linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    fig.suptitle(
        "Share of refinement's CheckM work that is removable",
        x=0.012,
        ha="left",
        fontsize=11.5,
        color=INK,
    )
    save(
        fig,
        "refinement_savings.png",
        NO_MW1
        + "  A bin is now keyed by a hash of its sorted contig ids, so a bin whose contents have been scored "
        "reuses that score - identical contigs necessarily give\nidentical completeness, so this is "
        "arithmetic rather than an approximation, and the cache is kept in work_files so --resume "
        "benefits. The size filter pays\nwhere bin sets arrive unfiltered: metaBAT2 emits its "
        "unbinnable remainder as a bin and on CAMI2 that reaches 159 Mb, so 0.3% of raw bins hold 11.8% "
        "of\nall the sequence. Inside bin_refinement it saves nothing, because a 50 kb - 20 Mb filter "
        "already runs on the inputs. An N50 floor was tried and dropped:\nthe worst N50 among passing "
        "bins is 2,025 bp, so a safe floor skips nothing and an effective one discards real genomes.",
    )


def fig_cost_landscape(agg) -> None:
    """What the pipeline actually spends, which is the context for every other figure here."""
    scores = read_json("legacy", "legacy_scores.json")
    if not scores:
        return
    by = scores["arms"]["modern"].get("by_step") or {}
    steps = [
        ("assembly", "assembly"),
        ("bin_refinement", "refinement"),
        ("binning", "binning"),
        ("read_qc", "read QC"),
    ]
    labels = [label for key, label in steps if by.get(key)]
    values = [by[key] for key, _l in steps if by.get(key)]
    share = [100 * v / sum(values) for v in values]

    fig, ax = plt.subplots(figsize=(8.6, 2.5))
    left = 0.0
    small = 0
    for i, (label, pct, value) in enumerate(zip(labels, share, values)):
        ax.barh([0], [pct], left=left, height=0.55, color=SERIES[i], linewidth=0)
        if pct > 4:
            ax.text(
                left + pct / 2,
                0,
                "%s\n%.1f cpu-h  (%.0f%%)" % (label, value, pct),
                ha="center",
                va="center",
                fontsize=8.5,
                color="white" if i == 0 else PAPER,
            )
        else:
            # Small segments get staggered labels with a leader, because two of them side by side at
            # one height print on top of each other.
            small += 1
            y = 0.30 if small % 2 else 0.46
            ax.annotate(
                "%s  %.1f cpu-h" % (label, value),
                (left + pct / 2, 0.28),
                xytext=(left + pct / 2 + 2.5, y),
                ha="left",
                fontsize=8,
                color=INK_SOFT,
                arrowprops={"arrowstyle": "-", "color": GRID, "linewidth": 0.9},
            )
        left += pct
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.5, 0.62)
    ax.set_yticks([])
    ax.set_xlabel("% of a default MetaWrap2 run's CPU time")
    ax.spines["left"].set_visible(False)
    fig.suptitle(
        "Share of CPU time by step, MetaWrap2 defaults",
        x=0.012,
        ha="left",
        fontsize=11.5,
        color=INK,
    )
    save(
        fig,
        "cost_landscape.png",
        "Three CAMI2 samples through MetaWrap2's defaults, 97.1 CPU-hours total. Binning is 1.8 of it, "
        "which is why the binner tiers in section 2 are affordable:\nmoving from the cheapest tier to "
        "the best measured one changes the whole run by about a third, not by the 20x the binning step "
        "itself changes by.",
    )


def main() -> int:
    style()
    os.makedirs(FIGURES, exist_ok=True)
    agg = load()
    fig_mw1(agg)
    fig_mw1_steps()
    fig_cost_landscape(agg)
    fig_binners(agg)
    fig_presets()
    fig_binners_alone(agg)
    fig_gpu(agg)
    fig_identity()
    fig_fidelity()
    fig_cohort()
    fig_quantifiers(agg)
    fig_refinement()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
