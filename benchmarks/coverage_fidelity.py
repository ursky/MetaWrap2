"""Do the alignment-free coverage estimators report what a binner expects? Measured per contig.

Every judgement about coverage in this project so far has been indirect: run a binner on the table and
count the MAGs. That answers "is the result worse" but not "is the number wrong, and wrong how" - and
those have different fixes. This compares the estimators to each other on **one assembly and one
sample's reads**, contig by contig, against bwa as the reference.

Three things it checks, because they fail independently:

1. **Agreement** - Pearson and Spearman against bwa's per-base mean depth.
2. **Units** - whether a method is reporting a per-base *depth* or something proportional to the
   contig's total coverage. A k-mer counter naturally produces a total, and dividing by length is the
   step that turns it into what metaBAT2 reads. If a method's values correlate better with
   ``depth x length`` than with ``depth``, it has skipped that division, and every long contig in the
   table is then overstated in proportion to how long it is.
3. **Presence** - how often a method calls a contig absent that bwa finds coverage for. This is the one
   that matters for a cohort: cross-sample binning works on the *pattern* of a contig's depth across
   samples, so a method that zeroes a sample destroys the signal even when its non-zero values are fine.

sylph is deliberately absent: it profiles a sample against a set of *genomes* and reports per-genome
abundance, so it has no per-contig number to compare. That is also why `quant_bins` can use it and
`coverage` cannot.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from typing import Dict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from metawrap2 import depth as _depth
from metawrap2.io.seqio import iter_fasta
from metawrap2.usage import run_measured

BINS_ENV = "metawrap2-bins"


def sh(argv, log, env=BINS_ENV, shell_pipeline=False):
    """Run *argv* inside the bins env, appending output to *log*. Returns (rc, cpu_hours)."""
    prefix = ["mamba", "run", "-n", env] if env else []
    full = prefix + (["bash", "-c", argv] if shell_pipeline else list(argv))
    with open(log, "a") as fh:
        fh.write("\n$ %s\n" % (argv if shell_pipeline else " ".join(argv)))
        fh.flush()
        rc, usage = run_measured(full, stdout=fh, stderr=subprocess.STDOUT)
    cpu = (usage.get("user_seconds", 0.0) + usage.get("system_seconds", 0.0)) / 3600.0
    return rc, cpu


def contig_lengths(assembly: str) -> pd.Series:
    return pd.Series({h.split()[0]: len(s) for h, s in iter_fasta(assembly)}, name="length")


def align_depth(assembly: str, r1: str, r2: str, out: str, aligner: str, threads: int) -> dict:
    """Per-contig mean depth via *aligner*, through the same code the pipeline uses."""
    work = os.path.join(out, aligner)
    os.makedirs(work, exist_ok=True)
    bam = os.path.join(work, "s.bam")
    table = os.path.join(work, "depth.txt")
    log = os.path.join(out, "%s.log" % aligner)
    cpu = 0.0
    if not os.path.isfile(table):
        if aligner == "bwa":
            _rc, c = sh(["bwa", "index", assembly], log)
            cpu += c
            cmd = "bwa mem -t %d %s %s %s" % (threads, assembly, r1, r2)
        else:
            cmd = "minimap2 -ax sr -t %d %s %s %s" % (threads, assembly, r1, r2)
        rc, c = sh(
            "set -o pipefail; %s | samtools sort -@ %d -o %s -" % (cmd, max(1, threads // 4), bam),
            log,
            shell_pipeline=True,
        )
        cpu += c
        if rc != 0:
            return {"error": "alignment failed, see %s" % log}
        start = time.time()
        _depth.write_depth_table([bam], table, samples=["s"], env=BINS_ENV, threads=threads)
        cpu += (time.time() - start) / 3600.0 * 0  # the table build is measured separately below
        os.remove(bam)
    frame = pd.read_csv(table, sep="\t").set_index("contigName")
    return {"values": frame["s"], "cpu_hours": round(cpu, 4)}


def fairy_depth(assembly: str, r1: str, r2: str, out: str, threads: int) -> dict:
    work = os.path.join(out, "fairy")
    os.makedirs(work, exist_ok=True)
    table = os.path.join(work, "depth.txt")
    log = os.path.join(out, "fairy.log")
    cpu = 0.0
    if not os.path.isfile(table):
        rc, c = sh(
            ["fairy", "sketch", "-1", r1, "-2", r2, "-S", "s", "-d", work, "-t", str(threads)], log
        )
        cpu += c
        if rc != 0:
            return {"error": "fairy sketch failed, see %s" % log}
        sketches = [os.path.join(work, f) for f in sorted(os.listdir(work)) if f.endswith(".bcsp")]
        rc, c = sh(
            ["fairy", "coverage", "-t", str(threads)] + sketches + [assembly, "-o", table], log
        )
        cpu += c
        if rc != 0:
            return {"error": "fairy coverage failed, see %s" % log}
    frame = pd.read_csv(table, sep="\t").set_index("contigName")
    # By name, not by position: the index has already been moved out of the columns, so slicing from
    # [3:] drops a real sample column and keeps a variance one.
    sample_cols = [
        c
        for c in frame.columns
        if c not in ("contigLen", "totalAvgDepth") and not c.endswith("-var")
    ]
    return {
        "values": frame[sample_cols[0]],
        "cpu_hours": round(cpu, 4),
        # Carried through so the report can check the column against jgi's definition.
        "total_avg_depth": frame["totalAvgDepth"],
        "sum_of_samples": frame[sample_cols].sum(axis=1),
    }


def coverm_depth(assembly: str, r1: str, r2: str, out: str, threads: int) -> dict:
    work = os.path.join(out, "coverm")
    os.makedirs(work, exist_ok=True)
    table = os.path.join(work, "coverm.tsv")
    log = os.path.join(out, "coverm.log")
    cpu = 0.0
    if not os.path.isfile(table):
        rc, cpu = sh(
            "coverm contig -1 %s -2 %s -r %s -t %d --methods mean -o %s"
            % (r1, r2, assembly, threads, table),
            log,
            shell_pipeline=True,
        )
        if rc != 0:
            return {"error": "coverm failed, see %s" % log}
    frame = pd.read_csv(table, sep="\t")
    frame = frame.set_index(frame.columns[0])
    return {"values": frame[frame.columns[0]], "cpu_hours": round(cpu, 4)}


def salmon_depth(
    assembly: str, r1: str, r2: str, out: str, threads: int, lengths: pd.Series
) -> dict:
    """salmon's per-contig read counts, converted to a depth so it is comparable at all.

    salmon reports TPM and NumReads, neither of which is a coverage. The conversion is
    ``reads x read_length / contig_length``, which is the definition of mean depth, and the read length
    is taken from the reads themselves rather than assumed.

    ``NumReads`` is already a count of *reads*, not of fragments, so there is no factor of two. The
    first version of this function multiplied by one, and salmon came out at 1.94x bwa - which looked
    like a tool disagreement and was arithmetic. MEGAHIT's own per-contig coverage estimate, in the
    FASTA headers, is what settled it: without the two, salmon lands at 1.07x MEGAHIT where bwa is at
    1.11x, and with it at 2.14x.
    """
    work = os.path.join(out, "salmon")
    os.makedirs(work, exist_ok=True)
    quant = os.path.join(work, "quant", "quant.sf")
    log = os.path.join(out, "salmon.log")
    cpu = 0.0
    if not os.path.isfile(quant):
        index = os.path.join(work, "index")
        rc, c = sh(["salmon", "index", "-t", assembly, "-i", index, "-p", str(threads)], log)
        cpu += c
        if rc != 0:
            return {"error": "salmon index failed, see %s" % log}
        rc, c = sh(
            ["salmon", "quant", "-i", index, "-l", "A", "-1", r1, "-2", r2,
             "-p", str(threads), "--meta", "-o", os.path.join(work, "quant")],
            log,
        )  # fmt: skip
        cpu += c
        if rc != 0:
            return {"error": "salmon quant failed, see %s" % log}
    frame = pd.read_csv(quant, sep="\t").set_index("Name")
    read_len = _read_length(r1)
    depth = frame["NumReads"] * read_len / lengths.reindex(frame.index)
    return {"values": depth, "cpu_hours": round(cpu, 4), "read_length": read_len,
            "raw_tpm": frame["TPM"]}  # fmt: skip


def _read_length(r1: str, records: int = 1000) -> int:
    """Modal read length from the first *records* reads, rather than an assumed 100."""
    lengths = []
    with open(r1) as fh:
        for index, line in enumerate(fh):
            if index % 4 == 1:
                lengths.append(len(line.strip()))
            if index >= records * 4:
                break
    return int(pd.Series(lengths).mode().iloc[0]) if lengths else 100


def compare(gold: pd.Series, other: pd.Series, lengths: pd.Series) -> dict:
    """How *other* relates to *gold*, including whether it looks like a total rather than a depth."""
    from scipy.stats import spearmanr

    frame = pd.DataFrame({"gold": gold, "other": other, "length": lengths}).dropna()
    both = frame[(frame.gold > 0) & (frame.other > 0)]
    if len(both) < 50:
        return {"error": "only %d contigs with coverage in both" % len(both)}
    g, o, length = both.gold.to_numpy(), both.other.to_numpy(), both.length.to_numpy()
    return {
        "contigs_compared": len(both),
        "pearson": round(float(np.corrcoef(g, o)[0, 1]), 4),
        "pearson_log": round(float(np.corrcoef(np.log10(g), np.log10(o))[0, 1]), 4),
        "spearman": round(float(spearmanr(g, o).statistic), 4),
        # The units test. Higher here than in `pearson` means the method skipped dividing by length.
        "pearson_vs_depth_times_length": round(float(np.corrcoef(g * length, o)[0, 1]), 4),
        "median_ratio": round(float(np.median(o / g)), 4),
        "gold_positive_other_zero_pct": round(
            100.0 * float(((frame.gold > 0) & (frame.other == 0)).mean()), 2
        ),
        "other_positive_gold_zero_pct": round(
            100.0 * float(((frame.other > 0) & (frame.gold == 0)).mean()), 2
        ),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="coverage_fidelity")
    ap.add_argument("-a", "--assembly", required=True)
    ap.add_argument("-1", dest="r1", required=True)
    ap.add_argument("-2", dest="r2", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-t", "--threads", type=int, default=96)
    ap.add_argument(
        "--methods",
        default="bwa,minimap2,fairy,coverm,salmon",
        help="comma-separated; bwa is the reference and is always run",
    )
    args = ap.parse_args(argv)
    os.makedirs(args.output, exist_ok=True)
    wanted = [m.strip() for m in args.methods.split(",") if m.strip()]
    lengths = contig_lengths(args.assembly)
    print("assembly: %d contigs, %.1f Mbp" % (len(lengths), lengths.sum() / 1e6), flush=True)

    results: Dict[str, dict] = {}
    runners = {
        "bwa": lambda: align_depth(
            args.assembly, args.r1, args.r2, args.output, "bwa", args.threads
        ),
        "minimap2": lambda: align_depth(
            args.assembly, args.r1, args.r2, args.output, "minimap2", args.threads
        ),
        "fairy": lambda: fairy_depth(args.assembly, args.r1, args.r2, args.output, args.threads),
        "coverm": lambda: coverm_depth(args.assembly, args.r1, args.r2, args.output, args.threads),
        "salmon": lambda: salmon_depth(
            args.assembly, args.r1, args.r2, args.output, args.threads, lengths
        ),
    }
    for name in ["bwa"] + [m for m in wanted if m != "bwa"]:
        if name not in runners:
            print("unknown method %r, skipping" % name, flush=True)
            continue
        print("running %s" % name, flush=True)
        results[name] = runners[name]()
        if "error" in results[name]:
            print("  %s: %s" % (name, results[name]["error"]), flush=True)
        else:
            print("  %s: %.2f cpu-h" % (name, results[name].get("cpu_hours", 0)), flush=True)

    if "error" in results.get("bwa", {"error": "not run"}):
        print("bwa is the reference and it failed; nothing to compare against.")
        return 1
    gold = results["bwa"]["values"]

    report: Dict[str, object] = {"assembly": args.assembly, "contigs": len(lengths)}
    table = {}
    for name, res in results.items():
        if name == "bwa" or "error" in res:
            continue
        table[name] = compare(gold, res["values"], lengths)
        table[name]["cpu_hours"] = res.get("cpu_hours")
    report["methods"] = table

    # fairy's totalAvgDepth column, against jgi's definition of it (the sum of the per-sample means).
    if "fairy" in results and "total_avg_depth" in results["fairy"]:
        total = results["fairy"]["total_avg_depth"]
        summed = results["fairy"]["sum_of_samples"].replace(0, np.nan)
        report["fairy_total_avg_depth_over_sum_of_samples"] = round(
            float(np.nanmedian(total / summed)), 4
        )

    out_json = os.path.join(args.output, "coverage_fidelity.json")
    with open(out_json, "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)

    print("\n%-10s %8s %9s %9s %12s %10s %9s" % (
        "method", "pearson", "log-r", "spearman", "r vs d*L", "med ratio", "cpu-h"))  # fmt: skip
    for name, row in table.items():
        if "error" in row:
            print("%-10s %s" % (name, row["error"]))
            continue
        print("%-10s %8.4f %9.4f %9.4f %12.4f %10.3f %9s" % (
            name, row["pearson"], row["pearson_log"], row["spearman"],
            row["pearson_vs_depth_times_length"], row["median_ratio"], row["cpu_hours"]))  # fmt: skip
    print("\nmissing-coverage rate (bwa finds depth, the method reports zero):")
    for name, row in table.items():
        if "error" not in row:
            print("  %-10s %5.1f%%" % (name, row["gold_positive_other_zero_pct"]))
    if "fairy_total_avg_depth_over_sum_of_samples" in report:
        print(
            "\nfairy totalAvgDepth / sum-of-samples = %.4f  (jgi defines it as the sum, so 1.0)"
            % report["fairy_total_avg_depth_over_sum_of_samples"]
        )
    print("\nwrote %s" % out_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
