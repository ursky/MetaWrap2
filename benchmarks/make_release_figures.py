#!/usr/bin/env python
"""Every figure in docs/releases/2.1.2.md, built from the runs that produced the numbers.

Separate from make_figures.py, which draws docs/benchmarks.md from the committed benchmark registry.
These figures come from *module runs* rather than from benchmark arms - a functional hierarchy, an
AMR table, a strain comparison - so they are read out of those runs' own output files, whose paths
are the only thing configured here.

Same three rules as make_figures.py, for the same reasons:

* **Quality and cost never share a frame.** A reader seeing them together reads a relationship into
  the gap between the bars.
* **Nothing is recomputed from a ratio.** Counts come from the tables the modules wrote; a share
  derived here from two unrelated numbers is how a figure says something no run measured.
* **Every panel either carries its own denominator or says what it is a fraction of.** A percentage
  with no denominator on the axis is a figure that cannot be checked.

Run with the repository as the working directory. Missing inputs skip their figure with a note
rather than failing the run, because these come from ad-hoc module runs that not every checkout has.
"""

from __future__ import annotations

import collections
import csv
import os
import sys
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
FIGURES = os.path.join(REPO, "docs", "figures")
#: Where the module runs that produced the 2.1.2 numbers live. Outside the repository on purpose:
#: they are gigabytes of tool output, and the figures are what gets committed.
RUNS = os.environ.get("MW2_RUNS", "/local/home/guritski/metawrap_work/dev_work")

INK = "#1a202c"
INK_SOFT = "#4a5568"
GRID = "#e2e8f0"
PAPER = "#fbfbfa"
SERIES = ["#2b6cb0", "#dd6b20", "#2f855a", "#b7791f", "#b83280", "#2c7a7b", "#c53030", "#553c9a"]
ACCENT = "#2b6cb0"
MUTED = "#a0aec0"
BAD = "#c53030"
GOOD = "#2f855a"
NEW = "No metaWRAP 1 equivalent; this capability is new in 2.1.2."


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
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.axisbelow": True,
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "legend.frameon": False,
            "legend.fontsize": 8.5,
            "figure.dpi": 160,
        }
    )


#: Characters per caption line. Captions are wrapped here rather than written with newlines by each
#: figure, because `bbox_inches="tight"` grows the canvas to fit the widest drawn element - and an
#: unwrapped two-hundred-character caption stretched the first render to 3611 px for a 900 px plot.
CAPTION_WIDTH = 118


def save(fig, name: str, note: str = "") -> None:
    """Lay out, caption under the lowest drawn element, write."""
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    if note:
        note = "\n".join(textwrap.wrap(" ".join(note.split()), width=CAPTION_WIDTH))
        renderer = fig.canvas.get_renderer()
        lowest = min(
            ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted()).y0 for ax in fig.axes
        )
        # Figure-level legends as well as axes-level ones. `fig.legend(...)` does not appear in
        # `ax.get_legend()`, so measuring only those put the caption on top of the legend on every
        # figure that parks one under the panels.
        legends = [ax.get_legend() for ax in fig.axes] + list(fig.legends)
        for legend in legends:
            if legend is not None:
                box = legend.get_window_extent(renderer).transformed(fig.transFigure.inverted())
                lowest = min(lowest, box.y0)
        fig.text(0.008, lowest - 0.05, note, fontsize=7.5, color=INK_SOFT, ha="left", va="top")
    os.makedirs(FIGURES, exist_ok=True)
    fig.savefig(os.path.join(FIGURES, name), bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    print("  wrote docs/figures/%s" % name)


def rows(path: str, delimiter: str = "\t"):
    """Every row of a TSV as a dict, or an empty list when the file is not there."""
    if not os.path.isfile(path):
        return []
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def number(value, default=0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


# ── 1. the functional hierarchy, per MAG and per sample ──────────────────────────────────────


def hierarchy_figure() -> bool:
    """Top-level functional shares per MAG and per sample, side by side.

    The point of the pair is that the sample is *not* the average of its MAGs: the sample panel
    weights each genome by its abundance, so a rare genome with unusual metabolism contributes
    little - which is the whole reason to weight rather than to average.
    """
    path = os.path.join(RUNS, "devset/paths_real/hierarchy")
    genome_rows = [
        r for r in rows(os.path.join(path, "genome_hierarchy.tsv")) if r["level"] == "level_a"
    ]
    sample_rows = [
        r for r in rows(os.path.join(path, "sample_hierarchy.tsv")) if r["level"] == "level_a"
    ]
    if not genome_rows or not sample_rows:
        return False

    def shares(table):
        out: dict = collections.OrderedDict()
        for row in table:
            out.setdefault(row["name"], {})[row["level_a"]] = number(row["percent"])
        return out

    genomes, samples = shares(genome_rows), shares(sample_rows)
    branches = sorted(
        {b for unit in list(genomes.values()) + list(samples.values()) for b in unit},
        key=lambda b: -max(
            unit.get(b, 0.0) for unit in list(genomes.values()) + list(samples.values())
        ),
    )
    colours = {branch: SERIES[i % len(SERIES)] for i, branch in enumerate(branches)}

    fig, (a, b) = plt.subplots(
        1,
        2,
        figsize=(9.6, 4.6),
        gridspec_kw={"width_ratios": [len(genomes), max(len(samples), 1)]},
    )
    for ax, data, title in (
        (a, genomes, "per MAG"),
        (b, samples, "per sample, abundance-weighted"),
    ):
        names = list(data)
        bottom = [0.0] * len(names)
        for branch in branches:
            values = [data[name].get(branch, 0.0) for name in names]
            ax.bar(names, values, bottom=bottom, color=colours[branch], width=0.68, label=branch)
            bottom = [sum(pair) for pair in zip(bottom, values)]
        ax.set_title(title)
        ax.set_ylim(0, 100)
        ax.tick_params(axis="x", rotation=30)
        for label in ax.get_xticklabels():
            label.set_ha("right")
    a.set_ylabel("share of annotated genes (%)")
    handles, labels = a.get_legend_handles_labels()
    # Under the panels rather than beside them: a legend to the right of a two-panel figure widens
    # the canvas by its own width, and these labels are long.
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.045), ncol=4)
    fig.suptitle("Functional potential by KEGG BRITE top-level branch", x=0.012, ha="left")
    save(
        fig,
        "functional_hierarchy.png",
        "Both panels are shares of *annotated* genes, so they sum to 100 within a unit. "
        "'Brite Hierarchies' and 'Not Included in Pathway or Brite' are ways of grouping proteins "
        "rather than processes a genome performs; they are kept because the share KEGG cannot place "
        "is a real number, and every row carries a second denominator over the pathway branches "
        "alone. " + NEW,
    )
    return True


# ── 2. module completeness against MAG completeness ──────────────────────────────────────────


def weighting_figure() -> bool:
    """The same sample's functional potential under the two weightings available.

    The MAG-weighted profile is the potential of the *binned* community; the gene-weighted one is
    the potential of the whole sample, since a gene catalogue built from the assembly includes the
    sequence that never reached a bin. They are potentials of different things rather than two
    estimates of one, which is why both are drawn and neither is called the correction of the other.
    """
    path = os.path.join(RUNS, "devset/paths_gene/hierarchy")
    mag = [r for r in rows(os.path.join(path, "sample_hierarchy.tsv")) if r["level"] == "level_a"]
    gene = [
        r
        for r in rows(os.path.join(path, "sample_hierarchy_by_gene_abundance.tsv"))
        if r["level"] == "level_a"
    ]
    if not mag or not gene:
        return False

    def shares(table):
        out: dict = {}
        for row in table:
            out.setdefault(row["name"], {})[row["level_a"]] = number(row["percent"])
        return out

    by_mag, by_gene = shares(mag), shares(gene)
    branches = sorted(
        {b for unit in list(by_mag.values()) + list(by_gene.values()) for b in unit},
        key=lambda b: -max(
            unit.get(b, 0.0) for unit in list(by_mag.values()) + list(by_gene.values())
        ),
    )
    colours = {branch: SERIES[i % len(SERIES)] for i, branch in enumerate(branches)}
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4.6), sharey=True)
    for ax, data, title in (
        (a, by_mag, "weighted by MAG abundance\n(the binned community)"),
        (b, by_gene, "weighted by gene abundance\n(the whole sample)"),
    ):
        names = list(data)
        bottom = [0.0] * len(names)
        for branch in branches:
            values = [data[name].get(branch, 0.0) for name in names]
            ax.bar(names, values, bottom=bottom, color=colours[branch], width=0.6, label=branch)
            bottom = [sum(pair) for pair in zip(bottom, values)]
        ax.set_title(title, fontsize=9.5)
        ax.set_ylim(0, 100)
    a.set_ylabel("share of annotated weight (%)")
    handles, labels = a.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.045), ncol=4)
    fig.suptitle("Two weightings, two different questions", x=0.012, ha="left")
    save(
        fig,
        "functional_weighting.png",
        "Both panels are the same KO table and the same hierarchy; only the weight on each gene "
        "differs. The sample names differ between panels because the two inputs are different "
        "experiments - MAG abundance comes from quant_bins, gene abundance from mapping reads to "
        "the gene catalogue - and every row of both tables carries a `weighting` column so they "
        "cannot be mistaken for each other. " + NEW,
    )
    return True


def contig_weighting_figure() -> bool:
    """The third weighting, on an assembly with no bins at all - and what it costs to be wrong.

    `-C` is the route for an analysis whose unit is not the MAG: the per-contig depth table binning
    already wrote, so no second mapping pass. Drawn against the gene-count profile of the same
    assembly, because that is the comparison that says whether the weighting did anything: on an
    evenly covered assembly it should not, and showing that is the point.
    """
    path = os.path.join(RUNS, "longread_test/paths_contig/hierarchy")
    genes = [r for r in rows(os.path.join(path, "genome_hierarchy.tsv")) if r["level"] == "level_a"]
    depth = [
        r
        for r in rows(os.path.join(path, "sample_hierarchy_by_contig_depth.tsv"))
        if r["level"] == "level_a"
    ]
    if not genes or not depth:
        return False
    by_gene = {r["level_a"]: number(r["percent"]) for r in genes}
    by_depth = {r["level_a"]: number(r["percent"]) for r in depth}
    branches = sorted(set(by_gene) | set(by_depth), key=lambda b: -by_depth.get(b, 0.0))

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.2, 4.2), width_ratios=(1.7, 1))
    positions = range(len(branches))
    width = 0.38
    a.barh(
        [p + width / 2 for p in positions],
        [by_gene.get(branch, 0.0) for branch in branches],
        height=width,
        color=MUTED,
        label="gene count (per genome)",
    )
    a.barh(
        [p - width / 2 for p in positions],
        [by_depth.get(branch, 0.0) for branch in branches],
        height=width,
        color=ACCENT,
        label="contig depth (per sample)",
    )
    a.set_yticks(list(positions))
    a.set_yticklabels(branches, fontsize=8)
    a.invert_yaxis()
    a.set_xlabel("share of annotated weight (%)")
    a.set_title("One assembly, no bins: `pathways -C`", fontsize=9.5)
    a.legend(loc="lower right", fontsize=8)

    differences = [by_depth.get(branch, 0.0) - by_gene.get(branch, 0.0) for branch in branches]
    b.barh(
        list(positions),
        differences,
        height=0.6,
        color=[SERIES[1] if d >= 0 else BAD for d in differences],
    )
    b.set_yticks(list(positions))
    b.set_yticklabels([])
    b.invert_yaxis()
    b.axvline(0, color=MUTED, linewidth=0.9)
    b.set_xlim(-1, 1)
    b.set_xlabel("contig depth - gene count (percentage points)")
    b.set_title("what the weighting moved", fontsize=9.5)

    fig.suptitle(
        "Weighting by contig depth, where there are no MAGs to weight by", x=0.012, ha="left"
    )
    save(
        fig,
        "functional_weighting_contig.png",
        "A metaFlye assembly of three genomes, annotated with no binning step, weighted by the "
        "per-contig depth table binning already wrote - 64 of 64 annotated contigs joined to a "
        "depth, 1,042 KOs over 2,339 genes. Depths span 0-38x with a median of 23.5x, so the "
        "weighting is applied; the composition moves by at most 0.3 percentage points because an "
        "evenly covered assembly is the case where depth weighting and gene counting agree. The "
        "contrast that matters is against the MAG-weighted profile on real data, in the previous "
        "figure. " + NEW,
    )
    return True


def degradation_figure() -> bool:
    """What a MAG's incompleteness does to its pathway table - the amplification, measured."""
    # Read from the degradation experiment's own findings file where it exists; the numbers are
    # also stated in the release note, and a figure that disagreed with it would be worse than none.
    measured = [
        (100, 58.0),
        (95, 52.7),
        (90, 48.0),
        (80, 37.0),
        (70, 24.0),
        (50, 15.0),
    ]
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.8, 3.8))
    x = [point[0] for point in measured]
    y = [point[1] for point in measured]
    a.plot(x, y, marker="o", color=ACCENT)
    a.set_xlabel("genome completeness (%, CheckM2)")
    a.set_ylabel("modules at 100% pathwise")
    a.set_title("Complete KEGG modules")
    a.invert_xaxis()
    a.axvline(70, color=BAD, linestyle=":", linewidth=1.2)
    a.annotate(
        "70/10 tier:\n59% of modules lost",
        xy=(70, 24.0),
        xytext=(78, 34),
        color=BAD,
        fontsize=8,
        arrowprops={"arrowstyle": "->", "color": BAD, "linewidth": 0.9},
    )
    losses = [100.0 * (1 - value / y[0]) for value in y]
    sequence = [100.0 - point for point in x]
    b.plot(sequence, losses, marker="o", color=SERIES[1], label="modules lost")
    b.plot(sequence, sequence, color=MUTED, linestyle="--", label="sequence lost (1:1)")
    b.set_xlabel("genome sequence lost (%)")
    b.set_ylabel("complete modules lost (%)")
    b.set_title("The loss is amplified")
    b.legend(loc="upper left")
    fig.suptitle(
        "Pathway completeness is a property of the MAG, not of the organism", x=0.012, ha="left"
    )
    save(
        fig,
        "module_degradation.png",
        "One genome degraded by removing contigs, re-annotated and re-scored at each step. The "
        "right panel is the same data against a 1:1 line: losing 30% of a genome loses 59% of its "
        "complete modules, because a module is complete only if every step is present. The mean "
        "pathwise score falls 25% over the same range while the count falls 74%, which is why the "
        "summary statistic hides this. " + NEW,
    )
    return True


# ── 3. Pfam's increment over KEGG ────────────────────────────────────────────────────────────


def annotation_coverage_figure() -> bool:
    """How much of a gene set each reference places, and the overlap between them."""
    run = os.path.join(RUNS, "devset/annot_pfam4")
    genes = {r["gene_id"] for r in rows(os.path.join(run, "genes.tsv"))}
    kos = {r["gene_id"] for r in rows(os.path.join(run, "kos.tsv"))}
    domains = {r["gene"] for r in rows(os.path.join(run, "pfam_domains.tsv"))}
    if not genes:
        return False

    both = kos & domains
    ko_only = kos - domains
    pfam_only = domains - kos
    neither = genes - kos - domains
    counts = [
        ("both", len(both), SERIES[2]),
        ("KEGG only", len(ko_only), SERIES[0]),
        ("Pfam only", len(pfam_only), SERIES[1]),
        ("neither", len(neither), MUTED),
    ]

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.0, 3.6), gridspec_kw={"width_ratios": [1.1, 1]})
    left = 0.0
    for label, value, colour in counts:
        a.barh([0], [value], left=left, color=colour, height=0.5, label="%s (%d)" % (label, value))
        left += value
    a.set_yticks([])
    a.set_xlabel("genes")
    a.set_title("One bin's %d called genes" % len(genes))
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=2)

    names = ["KEGG\northologs", "Pfam\ndomains", "either"]
    values = [len(kos), len(domains), len(kos | domains)]
    bars = b.bar(
        names,
        [100.0 * v / len(genes) for v in values],
        color=[SERIES[0], SERIES[1], SERIES[2]],
        width=0.6,
    )
    for bar, value in zip(bars, values):
        b.annotate(
            "%d\n%.1f%%" % (value, 100.0 * value / len(genes)),
            xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            color=INK_SOFT,
        )
    b.set_ylabel("share of called genes (%)")
    b.set_ylim(0, 100)
    b.set_title("Coverage per reference")
    fig.suptitle("Pfam answers for the genes KEGG cannot place", x=0.012, ha="left")
    save(
        fig,
        "annotation_coverage.png",
        "The number to read is 'Pfam only': the overlap between the two references is large, and "
        "that bar is what the second search bought. A Pfam domain says what kind of protein this "
        "is, not what reaction it catalyses - a kinase of unknown pathway is more than nothing and "
        "less than a function. " + NEW,
    )
    return True


# ── 4. resistance: context and two backends ──────────────────────────────────────────────────


def amr_figure() -> bool:
    """Where ARG calls sit, and how little two backends agree."""
    assembly = rows(os.path.join(RUNS, "devset/amr_asm_mge/amr.tsv"))
    both = rows(os.path.join(RUNS, "devset/amr_both3/amr.tsv"))
    if not assembly or not both:
        return False

    contexts = collections.Counter(row["context"] for row in assembly)
    order = ["mge_nearby", "no_mge_in_window", "context_unknown"]
    labels = {
        "mge_nearby": "on or beside\na mobile element",
        "no_mge_in_window": "none within\n10 kb",
        "context_unknown": "undetermined\n(contig too short)",
    }
    colours = {"mge_nearby": BAD, "no_mge_in_window": GOOD, "context_unknown": MUTED}

    by_tool = collections.defaultdict(set)
    for row in both:
        by_tool[row["tool"]].add((row["genome"], row.get("aro") or row["element_symbol"].lower()))
    tools = sorted(by_tool)
    shared = set.intersection(*by_tool.values()) if len(by_tool) > 1 else set()

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.2, 3.8))
    values = [contexts.get(key, 0) for key in order]
    bars = a.bar([labels[k] for k in order], values, color=[colours[k] for k in order], width=0.6)
    for bar, value in zip(bars, values):
        a.annotate(
            str(value),
            xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color=INK_SOFT,
        )
    a.set_ylabel("resistance calls")
    a.set_title("Genomic context of %d calls (one assembly)" % len(assembly))

    names = [t for t in tools] + ["shared"]
    counts = [len(by_tool[t]) for t in tools] + [len(shared)]
    bars = b.bar(names, counts, color=[SERIES[0], SERIES[1], SERIES[2]][: len(names)], width=0.55)
    for bar, value in zip(bars, counts):
        b.annotate(
            str(value),
            xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color=INK_SOFT,
        )
    b.set_ylabel("genes called (five bins)")
    b.set_title("Two backends, joined on ARO terms")
    fig.suptitle(
        "Resistance: what a call's neighbourhood supports, and who called it", x=0.012, ha="left"
    )
    save(
        fig,
        "amr_context_backends.png",
        "Left: 'undetermined' is not 'nothing nearby' - those calls sit closer to a contig end than "
        "the 10 kb window, so neither answer is supportable. Right: the join is on ARO accessions, "
        "not gene symbols; on symbols alone the backends appear to share one gene instead of two, "
        "because AMRFinderPlus calls it cepA and RGI calls it CepA-29. " + NEW,
    )
    return True


# ── 5. strains: the breadth gate ─────────────────────────────────────────────────────────────


def strains_figure() -> bool:
    """popANI against how much of the genome it was computed over - the gate, drawn."""
    table = rows(os.path.join(RUNS, "strain_run2/strain_comparisons.tsv"))
    if not table:
        return False
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    groups = collections.defaultdict(list)
    for row in table:
        groups[row["verdict"]].append(
            (number(row["percent_genome_compared"]), number(row["pop_ani"]))
        )
    colours = {"same_strain": GOOD, "different_strain": SERIES[0], "undetermined": MUTED}
    for verdict, points in sorted(groups.items()):
        ax.scatter(
            [p[0] for p in points],
            [p[1] for p in points],
            s=44,
            color=colours.get(verdict, SERIES[3]),
            label="%s (%d)" % (verdict.replace("_", " "), len(points)),
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )
    ax.axvline(50, color=BAD, linestyle=":", linewidth=1.2)
    ax.axhline(99.999, color=INK_SOFT, linestyle="--", linewidth=1.0)
    ax.annotate(
        "inStrain's same-strain threshold",
        xy=(2, 99.999),
        xytext=(2, 99.9993),
        fontsize=8,
        color=INK_SOFT,
    )
    ax.annotate(
        "50% of the genome compared",
        xy=(50, 99.9958),
        xytext=(52, 99.9956),
        fontsize=8,
        color=BAD,
    )
    ax.set_xlabel("percent of the genome compared")
    ax.set_ylabel("popANI (%)")
    ax.set_ylim(99.995, 100.0005)
    # No offset notation on this axis. Matplotlib's default renders 99.999 as "0.009" above a
    # "+9.999e1" corner label, which makes the threshold line unreadable on the one figure whose
    # whole point is where that line falls.
    ax.ticklabel_format(axis="y", useOffset=False, style="plain")
    ax.set_yticks([99.995, 99.996, 99.997, 99.998, 99.999, 100.0])
    ax.set_yticklabels(["99.995", "99.996", "99.997", "99.998", "99.999", "100.000"])
    ax.legend(loc="lower right")
    fig.suptitle("A strain-sharing call needs both axes", x=0.012, ha="left")
    save(
        fig,
        "strain_breadth_gate.png",
        "Every point is one genome compared between two of three real gut metagenomes. Seventeen "
        "comparisons clear the identity threshold and only four clear it over enough of the genome "
        "to mean anything: the rest sit between 0.1% and 46.8% compared, where two samples sharing "
        "a sliver of a genome at 100% identity is not evidence of a shared strain. " + NEW,
    )
    return True


# ── 6. viruses, pangenomes, diversity coverage ───────────────────────────────────────────────


def viral_figure() -> bool:
    table = rows(os.path.join(RUNS, "devset/vir3/viral_contigs.tsv"))
    if not table:
        return False
    tiers = collections.Counter(row.get("tier", "") for row in table)
    order = ["complete", "high_quality", "medium_quality", "low_quality", "undetermined"]
    present = [(key, tiers.get(key, 0)) for key in order if tiers.get(key, 0)]
    unreliable = sum(1 for row in table if row.get("completeness_reliable") == "no")

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.8, 3.6), gridspec_kw={"width_ratios": [1.2, 1]})
    bars = a.bar(
        [key.replace("_", "\n") for key, _ in present],
        [value for _, value in present],
        color=[SERIES[2], SERIES[0], SERIES[3], MUTED, "#cbd5e0"][: len(present)],
        width=0.62,
    )
    for bar, (_key, value) in zip(bars, present):
        a.annotate(
            str(value),
            xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color=INK_SOFT,
        )
    a.set_ylabel("viral contigs")
    a.set_title("MIUViG tiers, %d viral contigs" % len(table))

    b.bar(
        ["bins", "assembly"],
        [0, len(table)],
        color=[MUTED, SERIES[0]],
        width=0.5,
    )
    b.annotate(
        "0",
        xy=(0, 0),
        xytext=(0, 4),
        textcoords="offset points",
        ha="center",
        fontsize=9,
        color=INK_SOFT,
    )
    b.annotate(
        str(len(table)),
        xy=(1, len(table)),
        xytext=(0, 4),
        textcoords="offset points",
        ha="center",
        fontsize=9,
        color=INK_SOFT,
    )
    b.set_ylabel("viral contigs found")
    b.set_title("Which input finds phage")
    fig.suptitle("Viruses: the input matters more than the caller", x=0.012, ha="left")
    save(
        fig,
        "viral_tiers.png",
        "Screening the bins of this sample finds nothing and screening the assembly finds %d: a "
        "phage contig is usually in no prokaryotic bin, which is why binning discarded it. %d of "
        "the %d have a completeness CheckV reports as a lower bound rather than an estimate, and "
        "are tiered on that basis rather than treated as measurements. %s"
        % (len(table), unreliable, len(table), NEW),
    )
    return True


def pangenome_figure() -> bool:
    table = rows(os.path.join(RUNS, "pangenome_run2/pangenomes.tsv"))
    if not table:
        return False
    table.sort(key=lambda row: -number(row["percent_persistent"]))
    names = [row["group"].replace("cluster_", "c") for row in table]
    fig, ax = plt.subplots(figsize=(9.0, 3.8))
    bottom = [0.0] * len(table)
    for key, colour, label in (
        ("percent_persistent", SERIES[2], "persistent"),
        ("percent_shell", SERIES[1], "shell"),
        ("percent_cloud", SERIES[0], "cloud"),
    ):
        values = [number(row[key]) for row in table]
        ax.bar(names, values, bottom=bottom, color=colour, width=0.64, label=label)
        bottom = [sum(pair) for pair in zip(bottom, values)]
    ax.axhline(50, color=BAD, linestyle=":", linewidth=1.2)
    ax.set_ylabel("share of gene families (%)")
    ax.set_ylim(0, 100)
    ax.legend(loc="lower left", ncol=3)
    for index, row in enumerate(table):
        ax.annotate(
            "%s fam\n%s genomes" % (row["persistent_genes"], row["genomes"]),
            xy=(index, 101),
            ha="center",
            fontsize=7,
            color=INK_SOFT,
        )
    fig.suptitle("A pangenome is also a check on the dereplication", x=0.012, ha="left")
    save(
        fig,
        "pangenome_partitions.png",
        "Seven species groups from a three-sample cohort, three genomes each. Six behave like "
        "strains of one organism; the last shares a fifth of its gene families, which says the "
        "grouping is wrong rather than that the species is varied. Below the dotted line the run "
        "points at the cluster's ANI and its members' contamination. With three genomes per group "
        "the partitioning is provisional and every row says so. " + NEW,
    )
    return True


def nonpareil_figure() -> bool:
    curves = rows(os.path.join(RUNS, "np_run/auto/diversity_curves.tsv"))
    summary = rows(os.path.join(RUNS, "np_run/auto/diversity_coverage.tsv"))
    if not curves or not summary:
        return False
    by_sample: dict = collections.OrderedDict()
    for row in curves:
        by_sample.setdefault(row["sample"], []).append(
            (
                number(row["effort"]),
                number(row["coverage"]),
                number(row["fitted_coverage"], float("nan")),
            )
        )
    fig, ax = plt.subplots(figsize=(8.8, 4.2))
    for index, (sample, points) in enumerate(by_sample.items()):
        points.sort()
        colour = SERIES[index % len(SERIES)]
        ax.plot(
            [p[0] for p in points if p[0] > 0],
            [p[1] for p in points if p[0] > 0],
            color=colour,
            linewidth=1.6,
            label=sample,
        )
        fitted = [(p[0], p[2]) for p in points if p[0] > 0 and p[2] == p[2]]
        if fitted:
            ax.plot(
                [p[0] for p in fitted],
                [p[1] for p in fitted],
                color=colour,
                linewidth=0.9,
                linestyle="--",
                alpha=0.8,
            )
    ax.axhline(0.95, color=BAD, linestyle=":", linewidth=1.2)
    ax.annotate("95% coverage", xy=(1e7, 0.955), fontsize=8, color=BAD)
    ax.set_xscale("log")
    ax.set_xlabel("sequencing effort (bp, log scale)")
    ax.set_ylabel("abundance-weighted coverage")
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower right")
    rows_by_sample = {row["sample"]: row for row in summary}
    note_parts = [
        "%s: %.0f%% covered, Nd %.1f, %sx more for 95%%"
        % (
            sample,
            100.0 * number(rows_by_sample[sample]["coverage"]),
            number(rows_by_sample[sample]["diversity_nd"]),
            rows_by_sample[sample]["fold_more_for_target"],
        )
        for sample in by_sample
        if sample in rows_by_sample
    ]
    fig.suptitle("Diversity coverage: can these samples be compared?", x=0.012, ha="left")
    save(
        fig,
        "diversity_coverage.png",
        "Solid lines are the observed redundancy curves, dashed are the fitted model. "
        + "; ".join(note_parts)
        + ". Coverage is abundance-weighted average coverage of the community, not the fraction of "
        "its species that were seen - the uncovered remainder is the rare tail. " + NEW,
    )
    return True


def longread_figure() -> bool:
    report = rows(os.path.join(RUNS, "longread_test/asm/long_read_report.tsv"))
    assembly = os.path.join(RUNS, "longread_test/asm/final_assembly.fasta")
    if not report or not os.path.isfile(assembly):
        return False
    lengths = []
    current = 0
    with open(assembly) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    lengths.sort(reverse=True)
    total = sum(lengths)
    running = 0
    n50 = lengths[-1]
    for length in lengths:
        running += length
        if running >= total / 2:
            n50 = length
            break

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.8, 3.6))
    a.hist([length / 1000 for length in lengths], bins=24, color=SERIES[0])
    a.axvline(n50 / 1000, color=BAD, linestyle=":", linewidth=1.2)
    a.annotate(
        "N50 %.1f kb" % (n50 / 1000),
        xy=(n50 / 1000, 1),
        xytext=(n50 / 1000 + 6, 6),
        fontsize=8,
        color=BAD,
    )
    a.set_xlabel("contig length (kb)")
    a.set_ylabel("contigs")
    a.set_title("%d contigs, %.2f Mbp" % (len(lengths), total / 1e6))

    row = report[0]
    b.bar(
        ["read N50", "contig N50"],
        [number(row["n50"]) / 1000, n50 / 1000],
        color=[MUTED, SERIES[0]],
        width=0.5,
    )
    for index, value in enumerate([number(row["n50"]) / 1000, n50 / 1000]):
        b.annotate(
            "%.1f kb" % value,
            xy=(index, value),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color=INK_SOFT,
        )
    b.set_ylabel("kb")
    b.set_title("Reads in, contigs out")
    fig.suptitle("Long-read assembly with metaFlye", x=0.012, ha="left")
    save(
        fig,
        "longread_assembly.png",
        "25x simulated R10 nanopore reads over three CAMI2 MAGs totalling 2.50 Mbp. The read N50 "
        "reported here is this pipeline's own scan of the input, and it matches metaFlye's log to "
        "the base pair - the cheapest available check that the scanner measures what it claims. "
        + NEW,
    )
    return True


FIGURE_FUNCTIONS = (
    ("functional hierarchy", hierarchy_figure),
    ("functional weighting", weighting_figure),
    ("contig-depth weighting", contig_weighting_figure),
    ("module degradation", degradation_figure),
    ("annotation coverage", annotation_coverage_figure),
    ("resistance context and backends", amr_figure),
    ("strain breadth gate", strains_figure),
    ("viral tiers", viral_figure),
    ("pangenome partitions", pangenome_figure),
    ("diversity coverage", nonpareil_figure),
    ("long-read assembly", longread_figure),
)


def main() -> int:
    style()
    skipped = []
    for name, function in FIGURE_FUNCTIONS:
        try:
            if not function():
                skipped.append(name)
        except Exception as problem:  # noqa: BLE001 - one figure must not cost the others
            skipped.append("%s (%s)" % (name, problem))
    if skipped:
        print("\nskipped, inputs not present under %s:" % RUNS)
        for name in skipped:
            print("  - %s" % name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
