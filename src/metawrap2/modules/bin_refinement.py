"""MetaWrap2 bin_refinement module: consolidate several bin sets into one refined set.

Given any number of folders of bins produced from the same assembly (e.g. metaBAT2, SemiBin2,
COMEBin outputs), the default consolidator (:mod:`metawrap2.consolidation`) builds intersection, union and
difference candidates from overlapping bins, scores every candidate with CheckM2 (CheckM1 with
``--checkm1``), and keeps the best set of bins that share no contig.

``--mw1-refinement`` uses metaWRAP 1's method instead, the default up to 2.1.5: it hybridizes the
bin sets with Binning_refiner, scores every combined set, consolidates the best version of each bin
and dereplicates shared contigs. Its algorithms are in :mod:`metawrap2.refinement`.
``--consolidator dastool`` hands the input bin sets to DAS_Tool, which decides from single-copy
genes in one pass. Every route leaves one bin set in ``binsO``, and the remaining steps (re-scoring,
the -c/-x cut, plots, final layout, GUNC) are the same, so the routes can be compared on NC MAGs
recovered and cost. CheckM is in :mod:`metawrap2.checkm`.

External tool commands are in the COMMANDS block below and tunable numbers in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import itertools
import json
import os
import shlex
import shutil
from typing import Dict, List, Optional, Sequence, Tuple

from .. import checkm as _checkm
from .. import consolidation as _consolidation
from .. import gunc as _gunc
from .. import mimag as _mimag
from .. import recovery as _recovery
from .. import refinement as _refinement
from ..cliformat import bounded_float, make_parser, percent
from ..config import MODULE_ENVS, load_settings
from ..constants import (
    BIN_EXTENSION,
    CONTAMINATION_WEIGHT,
    CONTIGS_SUFFIX,
    FASTA_EXTENSIONS,
    MAX_BIN_SIZE,
    MIN_BIN_SIZE,
    STATS_SUFFIX,
    TSV_SUFFIX,
    bin_filename,
    is_refined_bins_dirname,
    is_unbinned_name,
    refined_bins_dirname,
)
from ..io.seqio import contig_id, fasta_files, iter_fasta
from ..progress import bar
from ..pyrun import PY
from ..utils import strip_trailing_slash
from ..validate import check_fasta
from ._common import (
    Raw,
    absolutize_paths,
    add_score_cache_arg,
    add_scorer_args,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    fill_command,
    finish_run,
    force,
    fs_replace_dir,
    listdir,
    make_checkpoint,
    memory_arg,
    open_score_cache,
    quoted_words,
    require_nonempty_dir,
    resolve_memory_gb,
    resolve_threads,
    run,
    scorer_notes,
    scorer_preflight,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["bin_refinement"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# MetaWrap2's own helpers run in the host interpreter (env=None), not the module's conda env,
# which holds CheckM and not metawrap2. See metawrap2.pyrun.
BINNING_REFINER = PY + " -m metawrap2.vendor.binning_refiner -1 {b1} -2 {b2} -o {out}"
BINNING_REFINER_3 = PY + " -m metawrap2.vendor.binning_refiner -1 {b1} -2 {b2} -3 {b3} -o {out}"
BINNING_REFINER_N = PY + " -m metawrap2.vendor.binning_refiner {bins} -o {out}"
PLOT_BINNING = PY + " -m metawrap2.scripts.plot_binning_results {comp} {cont} {stats}"
# --search_engine diamond: DAS_Tool's default, usearch, is not free software and is not in the env.
DAS_TOOL = (
    "DAS_Tool -i {tables} -l {labels} -c {assembly} -o {prefix} -t {threads}"
    " --search_engine diamond --score_threshold {score} --write_bins"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_MIN_COMPLETION = 70  # -c; bins must exceed this % completion to be "good"
DEFAULT_MAX_CONTAM = 10  # -x; bins must stay below this % contamination
# bin_refinement ignores bins outside this size range (bytes); see constants.py
MIN_SIZE = MIN_BIN_SIZE  # 50 kb
MAX_SIZE = MAX_BIN_SIZE  # 20 Mb
#: Upper limit on input bin sets for --mw1-refinement. It refines and scores combinations of the
#: sets, so cost grows quickly; this guards against accidental large runs and is not an algorithmic
#: limit. The default consolidator and dastool take any number of sets.
MAX_BIN_SETS = 6
MW1_MAX_BIN_SETS = MAX_BIN_SETS
#: Which bin-set combinations to build. "pairs-and-all" matches the original 2- and 3-set
#: behaviour (every pair plus the all-way combination); see refinement_plan().
COMBINATION_MODES = ("pairs-and-all", "pairs", "all")
DEFAULT_COMBINATIONS = "pairs-and-all"
#: How the bin sets are merged into one.
#:
#: "metawrap2" (default): intersection, union and difference candidates over cliques of
#: overlapping bins, all scored with CheckM, then disjoint bins selected by score with a lookahead
#: guard (metawrap2.consolidation). Unions let a bin be more complete than any input version of it.
#:
#: "metawrap1" (--mw1-refinement): metaWRAP 1's method, the default up to 2.1.5. Binning_refiner
#: intersections of every bin-set combination, each scored with CheckM, then the best version of
#: each bin chosen by pairwise 80% overlap, and shared contigs dereplicated.
#:
#: "dastool": score the input bins on single-copy genes and aggregate in one pass, a different
#: accuracy/cost trade-off. On per-sample bins, metaWRAP 1's method recovered 13 NC MAGs to DAS
#: Tool's 10, for 4.52 CPU-hours against 0.62.
CONSOLIDATORS = ("metawrap2", "dastool")
#: metaWRAP 1's method, selected with --mw1-refinement and not listed among --consolidator's choices.
MW1_CONSOLIDATOR = "metawrap1"
DEFAULT_CONSOLIDATOR = "metawrap2"
#: Older --consolidator values, still accepted: "hybrid" was the opt-in name of the default in
#: development builds, and "metawrap" named metaWRAP 1's method when it was the default.
CONSOLIDATOR_ALIASES = {"hybrid": "metawrap2", "metawrap": MW1_CONSOLIDATOR}
#: The default consolidator's candidate pool (see metawrap2.consolidation): graph edge at 10% shared
#: contigs, partners within 4x of each other's length, cliques of up to 3 bins, and derived candidates
#: between 200 kb and 20 Mb.
CANDIDATE_MIN_OVERLAP = _consolidation.MIN_OVERLAP_FRACTION
CANDIDATE_MAX_LENGTH_RATIO = _consolidation.MAX_LENGTH_RATIO
CANDIDATE_MAX_COMBO = _consolidation.MAX_COMBO
CANDIDATE_MIN_LENGTH = _consolidation.MIN_LENGTH
CANDIDATE_MAX_LENGTH = _consolidation.MAX_LENGTH
#: Directory for each kind of derived candidate, scored as one bin set each.
CANDIDATE_SET_PREFIX = "bins_candidates_"
#: Under work_files: the candidate pool, and which candidate became each consolidated bin.
POOL_FILE = "consolidation_pool.json"
CONSOLIDATED_BINS_FILE = "consolidated_bins.tsv"
#: DAS_Tool's default score threshold, a quality floor on the bins it keeps (like -c/-x).
DEFAULT_DASTOOL_SCORE = 0.5
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 bin_refinement",
        usage="metawrap2 bin_refinement [options] -o output_dir -b DIR1,DIR2[,DIR3...]",
        description="Consolidate any number of bin sets (from the same assembly) into a refined "
        "set. Contig names must be consistent across the bin folders.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    memory_arg(p)
    p.add_argument(
        "-c",
        "--completeness",
        type=percent(),
        default=DEFAULT_MIN_COMPLETION,
        help="minimum %% completion of bins (below 50 meets no MIMAG tier)",
    )
    p.add_argument(
        "-x",
        "--contamination",
        type=percent(),
        default=DEFAULT_MAX_CONTAM,
        help="maximum acceptable %% contamination of bins",
    )
    p.add_argument(
        "-b",
        "--bin-sets",
        action="append",
        default=[],
        dest="bin_sets",
        metavar="DIR1,DIR2,...",
        help="the bin sets to consolidate: folders of bins (.fa/.fasta) from the same assembly, "
        "comma-separated, or the option repeated once per folder. Order sets the internal names "
        "binsA, binsB, ... A value that names an existing folder is taken whole, so a folder "
        "whose name contains a comma can be given with -b on its own",
    )
    # Older spellings: --bins DIR (repeated) and --bin_sets mean -b. -A..-F named up to six sets
    # positionally; they still work but cannot be mixed with -b.
    for alias in ("--bin_sets", "--bins"):
        p.add_argument(alias, action="append", dest="bin_sets", help=argparse.SUPPRESS)
    p.add_argument("-A", "--bins-a", help=argparse.SUPPRESS)
    for letter in "BCDEF":
        p.add_argument("-" + letter, "--bins-" + letter.lower(), help=argparse.SUPPRESS)
    p.add_argument(
        "--refine-combinations",
        choices=COMBINATION_MODES,
        default=None,
        help="with --mw1-refinement, which bin-set combinations to build with Binning_refiner: "
        "'pairs-and-all' (default) every pair plus the all-way combination; "
        "'pairs' every pair only (cheaper); 'all' every subset of size >=2 "
        "(the number grows exponentially). Has no effect with 2 or 3 bin sets. Implies "
        "--mw1-refinement unless --consolidator is given",
    )
    p.add_argument(
        "--skip-refinement",
        dest="refine",
        action="store_false",
        help="don't use Binning_refiner to make refined bins from binner combinations. Implies "
        "--mw1-refinement unless --consolidator is given",
    )
    p.add_argument(
        "--score-everything",
        action="store_true",
        help="score every bin with CheckM, including ones excluded by size or fragmentation. "
        "Slower; use it to check that the filter is correct",
    )
    p.add_argument(
        "--skip-checkm",
        dest="run_checkm",
        action="store_false",
        help="don't run CheckM to assess bins. Implies --mw1-refinement unless --consolidator "
        "is given",
    )
    p.add_argument(
        "--skip-gunc",
        dest="run_gunc",
        action="store_false",
        help="don't run GUNC on the final bins. This skips the chimerism check, which CheckM "
        "does not do, and the need for the GUNC database",
    )
    # `score` spells it --no-gunc; accept that here too.
    p.add_argument("--no-gunc", dest="run_gunc", action="store_false", help=argparse.SUPPRESS)
    p.add_argument(
        "--skip-consolidation",
        dest="consolidate",
        action="store_false",
        help="don't choose the best version of each bin; just pick the best folder. Implies "
        "--mw1-refinement unless --consolidator is given",
    )
    amb = p.add_mutually_exclusive_group()
    amb.add_argument(
        "--keep-ambiguous",
        dest="dereplicate",
        action="store_const",
        const="keep",
        help="for contigs in more than one bin, keep them in all bins (--mw1-refinement; "
        "implies it unless --consolidator is given)",
    )
    amb.add_argument(
        "--remove-ambiguous",
        dest="dereplicate",
        action="store_const",
        const="remove",
        help="for contigs in more than one bin, remove them from all bins (--mw1-refinement; "
        "implies it unless --consolidator is given)",
    )
    p.set_defaults(dereplicate="partial")  # default: keep only in the best bin
    p.add_argument(
        "--checkm-jobs",
        type=int,
        default=0,
        metavar="N",
        help="number of bin sets to score with CheckM in parallel (default 0: as many as memory "
        "allows; 1: one at a time)",
    )
    p.add_argument(
        "--checkm2-ttable",
        type=int,
        choices=(4, 11),
        default=None,
        metavar="N",
        help="run CheckM2 under this translation table only (4 or 11) instead of trying both. "
        "Ignored with --checkm1",
    )
    add_scorer_args(p)
    add_score_cache_arg(p)
    p.add_argument(
        "--consolidator",
        default=None,
        metavar="{%s}" % ",".join(CONSOLIDATORS),
        help="how to merge the bin sets into one set: 'metawrap2' (default) adds unions and "
        "differences of overlapping bins to their intersections, scores every candidate with "
        "CheckM and keeps the best non-overlapping set; 'dastool' scores the input bins on "
        "single-copy genes and aggregates in one pass, and needs -a and the metawrap2-binqc env "
        "(created by `metawrap2 install-env bin_refinement`)",
    )
    p.add_argument(
        "--mw1-refinement",
        action="store_true",
        help="consolidate with metaWRAP 1's method, the default up to 2.1.5: Binning_refiner "
        "intersections of every bin-set combination, the best version of each bin chosen by "
        "80%% overlap, and shared contigs dereplicated",
    )
    p.add_argument(
        "--contamination-weight",
        type=bounded_float(0, None),
        default=None,
        metavar="W",
        help="rank candidate bins by completeness - W * contamination (default %g). A lower W "
        "favours more complete bins over cleaner ones. Refused with --mw1-refinement and "
        "--consolidator dastool, which rank with the fixed weight %g"
        % (CONTAMINATION_WEIGHT, CONTAMINATION_WEIGHT),
    )
    p.add_argument(
        "-a",
        "--assembly",
        help="assembly the bins came from (required by --consolidator dastool, which reads the "
        "contigs itself; the default and --mw1-refinement do not need it)",
    )
    p.add_argument(
        "--dastool-score",
        type=bounded_float(0, 1),
        default=DEFAULT_DASTOOL_SCORE,
        help="DAS_Tool --score_threshold: minimum bin score to keep",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    args = p.parse_args(argv)
    if args.consolidator is not None and args.consolidator not in (
        CONSOLIDATORS + (MW1_CONSOLIDATOR,) + tuple(CONSOLIDATOR_ALIASES)
    ):
        p.error(
            "argument --consolidator: invalid choice: %r (choose from %s)"
            % (args.consolidator, ", ".join(repr(c) for c in CONSOLIDATORS))
        )
    return args


def _mw1_only_options(args) -> List[str]:
    """The options given in *args* that only metaWRAP 1's method uses."""
    found = []
    if not args.run_checkm:
        found.append("--skip-checkm")
    if not args.refine:
        found.append("--skip-refinement")
    if not args.consolidate:
        found.append("--skip-consolidation")
    if args.refine_combinations is not None:
        found.append("--refine-combinations")
    if args.dereplicate == "keep":
        found.append("--keep-ambiguous")
    elif args.dereplicate == "remove":
        found.append("--remove-ambiguous")
    return found


def resolve_consolidator(args) -> str:
    """Set ``args.consolidator`` to the method this run uses, and return it.

    The result is ``metawrap2`` (the default), ``metawrap1`` (metaWRAP 1's method) or ``dastool``.

    * An old name in --consolidator is mapped (``hybrid`` to metawrap2, ``metawrap`` to metaWRAP 1's
      method), with a note.
    * --mw1-refinement selects metaWRAP 1's method. It conflicts with --consolidator metawrap2 or
      dastool.
    * Without --consolidator, an option that only metaWRAP 1's method uses (--skip-checkm,
      --skip-refinement, --skip-consolidation, --refine-combinations, --keep-ambiguous,
      --remove-ambiguous) selects it, with a note, so those options work as they did up to 2.1.5.
      With an explicit --consolidator metawrap2 they are an error.
    * --contamination-weight is an error with any method other than metawrap2.
    """
    given = args.consolidator
    if given in CONSOLIDATOR_ALIASES:
        mapped = CONSOLIDATOR_ALIASES[given]
        comm(
            "--consolidator %s is an older name; using %s."
            % (
                given,
                "--mw1-refinement" if mapped == MW1_CONSOLIDATOR else "--consolidator " + mapped,
            )
        )
        given = mapped
    if args.mw1_refinement:
        if given in CONSOLIDATORS:
            error(
                "--mw1-refinement and --consolidator %s cannot be used together: each picks a "
                "consolidator.\n"
                "Drop --consolidator to use metaWRAP 1's method, or drop --mw1-refinement to use "
                "%s." % (given, given)
            )
        given = MW1_CONSOLIDATOR
    mw1_options = _mw1_only_options(args)
    if given is None:
        if mw1_options:
            comm(
                "%s applies to metaWRAP 1's refinement, so this run uses it (--mw1-refinement)."
                % " and ".join(mw1_options)
            )
            given = MW1_CONSOLIDATOR
        else:
            given = DEFAULT_CONSOLIDATOR
    elif given == "metawrap2" and mw1_options:
        error(
            "--consolidator metawrap2 and %s cannot be used together: %s applies to metaWRAP 1's "
            "refinement only.\n"
            "Drop --consolidator metawrap2 to use metaWRAP 1's method, or drop %s."
            % (mw1_options[0], mw1_options[0], mw1_options[0])
        )
    if args.contamination_weight is not None and given != "metawrap2":
        label = "--mw1-refinement" if given == MW1_CONSOLIDATOR else "--consolidator " + given
        error(
            "--contamination-weight cannot be used with %s, which ranks bins with the fixed "
            "weight %g.\n"
            "Drop --contamination-weight, or drop %s to use the default consolidator."
            % (label, CONTAMINATION_WEIGHT, label)
        )
    if args.contamination_weight is None:
        args.contamination_weight = CONTAMINATION_WEIGHT
    if args.refine_combinations is None:
        args.refine_combinations = DEFAULT_COMBINATIONS
    args.consolidator = given
    return given


def _split_bin_sets(values: Sequence[str]) -> List[str]:
    """The folders named by -b values: each comma-separated, unless the value is a folder."""
    found: List[str] = []
    for value in values:
        if os.path.isdir(os.path.expanduser(value)) or "," not in value:
            found.append(value)
        else:
            found.extend(part.strip() for part in value.split(",") if part.strip())
    return found


def collect_bin_sets(args) -> List[str]:
    """The input bin-set directories, in order: from -b, or from the older -A..-F.

    Order sets each label (binsA, binsB, ...), and therefore the refinement plan and the bin
    numbering in each refined set, so it is part of the interface.
    """
    lettered = [
        getattr(args, "bins_" + letter) for letter in "abcdef" if getattr(args, "bins_" + letter)
    ]
    listed = _split_bin_sets(getattr(args, "bin_sets", None) or [])
    if lettered and listed:
        error(
            "Give the bin sets one way: -b DIR1,DIR2,... or the older -A/-B/-C..., not both.\n"
            "For example: -b %s" % ",".join(lettered + listed)
        )
    ordered = [os.path.abspath(os.path.expanduser(p)) for p in (lettered or listed)]
    if not ordered:
        error("No bin sets given. Pass them with -b DIR1,DIR2,... (or -b once per folder).")
    seen: Dict[str, str] = {}
    for path in ordered:
        key = os.path.realpath(path)
        if key in seen:
            error(
                "The bin set %s is given twice%s. Give each folder once."
                % (path, "" if seen[key] == path else " (also as %s)" % seen[key])
            )
        seen[key] = path
    if not dry_run():
        # A dry run plans a study whose bin folders do not exist yet.
        missing = [path for path in ordered if not os.path.isdir(path)]
        if missing:
            error(
                "Not a folder: %s.\nEach bin set is a folder of bins. Paths given to -b are "
                "separated by commas, so a path with a comma in it must name an existing folder "
                "and be given with -b on its own." % ", ".join(missing)
            )
    return ordered


def drop_empty_bin_sets(bin_sets: Sequence[str]) -> List[str]:
    """The bin sets that hold at least one bin, warning about each one that holds none.

    A binner with nothing to cluster writes an empty folder and exits 0 (metaBAT2 on a small
    assembly does), and the other sets can still be refined. Only an existing folder with no FASTA
    is dropped; a missing path is left for the input checks. If every set is empty this returns an
    empty list and the caller writes an empty final set.
    """
    kept: List[str] = []
    for index, path in enumerate(bin_sets):
        if os.path.isdir(path) and not any(
            name.endswith(tuple(FASTA_EXTENSIONS)) for name in os.listdir(path)
        ):
            warning("bin set %d (%s) has no bins in it; refining without it." % (index + 1, path))
            continue
        kept.append(path)
    if not kept:
        warning(
            "Every bin set given is empty (%s); writing an empty final bin set.\n"
            "Check the binning logs. A binner that makes no bins usually had too few contigs "
            "above its length cut-off, or no coverage for them." % ", ".join(bin_sets)
        )
    elif len(kept) < len(bin_sets):
        comm(
            "Refining %d of %d bin sets; they are relabelled %s in order."
            % (len(kept), len(bin_sets), ", ".join(bin_set_label(i) for i in range(len(kept))))
        )
    return kept


def write_empty_result(out: str, completeness: float, contamination: float) -> str:
    """Write the final bin set for a run with no bins: an empty folder and a header-only table."""
    final = os.path.join(out, refined_bins_dirname(completeness, contamination))
    os.makedirs(final, exist_ok=True)
    with open(final + STATS_SUFFIX, "w") as fh:
        fh.write("\t".join(_checkm.STATS_HEADER) + "\n")
    return final


def bin_set_label(index: int) -> str:
    """The internal name of the *index*-th input bin set: binsA ... binsZ, then binsAA, binsAB, ..."""
    letters = ""
    index += 1
    while index:
        index, rest = divmod(index - 1, 26)
        letters = chr(ord("A") + rest) + letters
    return "bins" + letters


def combined_label(labels: Sequence[str]) -> str:
    """The name of the set refined from *labels*, e.g. (binsA, binsC) -> binsAC."""
    return "bins" + "".join(label[len("bins") :] for label in labels)


#: The original refinement plan for 2 and 3 inputs, kept unchanged.
#:
#: The build order and the orientation of each pair both matter. Binning_refiner numbers output
#: bins by input folder order, so binsBC built as (C, B), as this module has always done, numbers
#: bins differently from (B, C). Changing this changes results for existing 2- and 3-binner runs.
_LEGACY_PLANS: Dict[int, List[Tuple[Tuple[str, ...], str]]] = {
    2: [(("binsA", "binsB"), "binsAB")],
    3: [
        (("binsA", "binsB"), "binsAB"),
        (("binsC", "binsB"), "binsBC"),
        (("binsA", "binsC"), "binsAC"),
        (("binsA", "binsB", "binsC"), "binsABC"),
    ],
}


def refinement_plan(
    labels: Sequence[str], combinations: str = DEFAULT_COMBINATIONS
) -> List[Tuple[Tuple[str, ...], str]]:
    """The bin-set combinations to refine, as [(input labels, output label), ...].

    For two or three inputs this returns :data:`_LEGACY_PLANS` unchanged, so those runs match
    earlier versions. For more inputs:

    ``pairs-and-all`` (default)
        every pair, then the combination of all of them, as in the 3-input case (3 pairs plus the
        triple).
    ``pairs``
        every pair only. Cheaper: the number of CheckM runs grows quadratically, not
        exponentially.
    ``all``
        every subset of size >= 2. Thorough and expensive: 2**n - n - 1 combinations.
    """
    if combinations not in COMBINATION_MODES:
        raise ValueError(
            "combinations must be one of %s, got %r" % (", ".join(COMBINATION_MODES), combinations)
        )
    labels = list(labels)
    if len(labels) < 2:
        return []
    if combinations == DEFAULT_COMBINATIONS and len(labels) in _LEGACY_PLANS:
        # Only valid for the canonical A/B/C labels, which is how main() names them.
        legacy = _LEGACY_PLANS[len(labels)]
        if labels == [bin_set_label(i) for i in range(len(labels))]:
            return list(legacy)

    plan: List[Tuple[Tuple[str, ...], str]] = []
    if combinations == "all":
        sizes = range(2, len(labels) + 1)
    else:
        sizes = range(2, 3)
    for size in sizes:
        for combo in itertools.combinations(labels, size):
            plan.append((combo, combined_label(combo)))
    if combinations == DEFAULT_COMBINATIONS and len(labels) > 2:
        all_way = tuple(labels)
        plan.append((all_way, combined_label(all_way)))
    return plan


def _fix_contig_naming(path: str) -> bool:
    """Normalise a bin fasta's headers in place. Returns True if the file was rewritten.

    The header is reduced to the contig id (:func:`contig_id`), dropping description text that
    binners append (current metaBAT2 writes ``total_depth=.. sample_depths=..``, which made its bins
    match nothing). Any remaining '=' becomes '_', as in fix_config_naming.py, because '=' breaks
    some downstream tools.

    Files with clean headers are not rewritten, which avoids needless I/O and keeps their mtimes.
    """
    records = list(iter_fasta(path))
    fixed = [(contig_id(h).replace("=", "_"), seq) for h, seq in records]
    if all(new == old for (new, _), (old, _) in zip(fixed, records)):
        return False
    with open(path, "w") as out:
        out.writelines(">%s\n%s\n" % (header, seq) for header, seq in fixed)
    return True


def _copy_bins(src: str, dest: str) -> int:
    """Copy size-filtered bins from *src* into *dest* as bin-name.fa. Returns count kept.

    The unbinned catch-all is skipped by name. `binning` writes CONCOCT's leftovers to
    ``concoct_bins/unbinned.fasta`` so users can see what no binner claimed, and the AMBER scorer
    also skips it. Copied in, it would go through CheckM, every Binning_refiner combination and
    dereplication as if it were a genome. On CAMI2 GI the size filter happened to drop it (69.8 MB,
    over 20 MB), but a smaller catch-all would pass.

    What the size filter drops is summarised separately for each end. Bins over the ceiling hold
    too much sequence to be one genome. Bins under the floor are fragments that might have grown
    with a different binner in the set. On a CAMI2 CONCOCT set the filter removed 10 of 99 bins
    holding 41% of the sequence (7 over the ceiling, the largest 131.8 Mb), which accounts for most
    of the recall refinement appears to lose. Percentages use file size, as the filter does.
    """
    ensure_dir(dest)
    n = 0
    dropped: dict[str, list[int]] = {"over": [], "under": []}
    total = 0
    for name in sorted(os.listdir(src)):
        full = os.path.join(src, name)
        if not os.path.isfile(full):
            continue
        if is_unbinned_name(name):
            comm("Skipping %s: it holds the unbinned contigs, not a bin" % full)
            continue
        size = os.path.getsize(full)  # portable (no `stat -c`)
        total += size
        if MIN_SIZE < size < MAX_SIZE:
            base = os.path.splitext(name)[0]
            shutil.copy(full, os.path.join(dest, base + BIN_EXTENSION))
            n += 1
        else:
            dropped["over" if size >= MAX_SIZE else "under"].append(size)
            comm("Skipping %s because the bin size is not between 50kb and 20Mb" % full)
    for end, label in (("over", "above the 20Mb ceiling"), ("under", "below the 50kb floor")):
        sizes = dropped[end]
        if sizes:
            comm(
                "  %d bin(s) %s, holding %.1f%% of the set's sequence (largest %.1fMb)"
                % (len(sizes), label, 100.0 * sum(sizes) / total, max(sizes) / 1e6)
            )
    return n


def _count_bins(path: str) -> int:
    return (
        len([f for f in os.listdir(path) if f.endswith(BIN_EXTENSION)])
        if os.path.isdir(path)
        else 0
    )


def write_contig_table(bins_dir: str, path: str) -> int:
    """Write *bins_dir*'s contig-to-bin membership to *path* as TSV. Returns the row count.

    Used for the table beside each final bin set (read by Anvio and similar tools) and for
    DAS_Tool's per-bin-set input. One writer keeps the format identical; DAS_Tool given a
    differently formatted table consolidates nothing without an error.
    """
    rows = 0
    with open(path, "w") as out_fh:
        for name in sorted(os.listdir(bins_dir)):
            if not name.endswith(BIN_EXTENSION):
                continue
            bin_name = os.path.splitext(name)[0]
            for header, _seq in iter_fasta(os.path.join(bins_dir, name)):
                out_fh.write("%s\t%s\n" % (contig_id(header), bin_name))
                rows += 1
    return rows


def assembly_for_bins(assembly: str, work: str) -> str:
    """The assembly to hand a consolidator, with headers that match the copied bins.

    The bins have been through :func:`_fix_contig_naming` (header reduced to the contig id, '='
    replaced by '_'). A consolidator matches assembly contig names against the bin tables, so
    unnormalised names give an empty consolidation without an error.

    Records with no sequence are also dropped. DAS_Tool reads the assembly with R's ``fread``; an
    empty record gives a one-field line, the read stops, and every later contig is reported "not
    found in assembly" (36,854 of them). CAMI2's gold-standard assembly has two such records;
    `binning` already removes them from its own working copy.

    A copy is written only if something changed (not needed for metaSPAdes or MEGAHIT names with
    no empty records), since assemblies can be gigabytes.
    """
    renamed = [(contig_id(h).replace("=", "_"), h, seq) for h, seq in iter_fasta(assembly)]
    empty = [new for new, _old, seq in renamed if not seq]
    if not empty and all(new == contig_id(old) for new, old, _seq in renamed):
        return assembly
    ensure_dir(work)
    normalized = os.path.join(work, "assembly_normalized.fa")
    reasons = []
    if any(new != contig_id(old) for new, old, _seq in renamed):
        reasons.append("headers need the same normalisation the bins got")
    if empty:
        reasons.append(
            "%d record(s) have no sequence, which DAS_Tool's assembly reader stops on" % len(empty)
        )
    comm("writing %s: %s" % (os.path.basename(normalized), "; ".join(reasons)))
    with open(normalized, "w") as fh:
        fh.writelines(">%s\n%s\n" % (new, seq) for new, _old, seq in renamed if seq)
    return normalized


def _run_dastool(args, env: Optional[str], out: str, labels: Sequence[str]) -> None:
    """Consolidate the bin sets in *out* with DAS_Tool, leaving the result as ``binsO``.

    binsO is where the built-in consolidation also writes, so the later steps (CheckM re-scoring,
    the -c/-x cut, ranking plots, final layout, GUNC) are the same for both consolidators.
    """
    work = os.path.join(out, "work_files")
    ensure_dir(work)
    tables = []
    for label in labels:
        table = os.path.join(work, label + ".dastool.tsv")
        if not dry_run():
            rows = write_contig_table(os.path.join(out, label), table)
            comm("%s contributes %d contigs to DAS_Tool" % (label, rows))
        tables.append(table)

    assembly = args.assembly if dry_run() else assembly_for_bins(args.assembly, work)
    prefix = os.path.join(work, "dastool", "DASTool")
    ensure_dir(os.path.dirname(prefix))
    run(
        fill_command(
            DAS_TOOL,
            tables=",".join(tables),
            labels=",".join(labels),
            assembly=assembly,
            prefix=prefix,
            threads=args.threads,
            score=args.dastool_score,
        ),
        env=env,
        tool="DAS_Tool",
    )
    if dry_run():
        return
    produced = prefix + "_DASTool_bins"
    if not os.path.isdir(produced):
        error(
            "DAS_Tool did not produce %s.\n"
            "Its log is in %s. This happens when DAS_Tool finds no single-copy genes, usually "
            "because the contig names in the bins do not match the assembly given "
            "with -a." % (os.path.basename(produced), os.path.dirname(prefix))
        )
    binsO = os.path.join(out, "binsO")
    ensure_dir(binsO)
    # Renamed to bin_001.fasta, matching every MetaWrap2 bin set and the .stats.tsv keys. DAS_Tool's
    # names record each bin's input set, so the mapping is written to dastool_bin_names.tsv.
    kept = 0
    provenance = []
    for name in sorted(os.listdir(produced)):
        if os.path.splitext(name)[1] in FASTA_EXTENSIONS:
            kept += 1
            renamed = bin_filename(kept)
            shutil.copy(os.path.join(produced, name), os.path.join(binsO, renamed))
            provenance.append((os.path.splitext(renamed)[0], name))
    with open(os.path.join(work, "dastool_bin_names.tsv"), "w") as fh:
        fh.write("bin\tdastool_name\n")
        fh.writelines("%s\t%s\n" % row for row in provenance)
    comm("DAS_Tool kept %d bins (score threshold %s)" % (kept, args.dastool_score))
    if kept == 0:
        error(
            "DAS_Tool kept no bins at --dastool-score %s. Its per-bin scores are in "
            "%s_DASTool_summary.tsv; lower the threshold if the best score there is acceptable."
            % (args.dastool_score, prefix)
        )


def _candidate_set_dirs(out: str) -> List[str]:
    """The derived-candidate sets of a default-consolidator run present in *out*, sorted."""
    return [
        name
        for name in sorted(listdir(out))
        if name.startswith(CANDIDATE_SET_PREFIX) and os.path.isdir(os.path.join(out, name))
    ]


def _build_candidate_pool(out: str, labels: Sequence[str]) -> None:
    """Build the default consolidator's candidate pool from the input bin sets in *out*.

    The input bins are candidates where they are, in binsA, binsB, ... Each kind of derived
    candidate is written as its own bin set (``bins_candidates_union`` and so on) so the CheckM step
    scores it with the input sets, in parallel under ``--checkm-jobs``. Contig sequences come from
    the input bins, so no assembly is needed. The pool, with each candidate's contigs, is recorded
    in ``work_files/consolidation_pool.json`` for the consolidation step.
    """
    bins: Dict[str, frozenset] = {}
    files: Dict[str, str] = {}
    sequences: Dict[str, str] = {}
    for label in labels:
        for path in fasta_files(os.path.join(out, label)):
            name = "%s/%s" % (label, os.path.splitext(os.path.basename(path))[0])
            contigs = []
            for header, seq in iter_fasta(path):
                contig = contig_id(header)
                contigs.append(contig)
                sequences.setdefault(contig, seq)
            if contigs:
                bins[name] = frozenset(contigs)
                files[name] = os.path.relpath(path, out)
    lengths = {contig: len(seq) for contig, seq in sequences.items()}
    pool = _consolidation.build_candidates(
        bins,
        lengths,
        max_combo=CANDIDATE_MAX_COMBO,
        min_overlap=CANDIDATE_MIN_OVERLAP,
        max_length_ratio=CANDIDATE_MAX_LENGTH_RATIO,
        min_length=CANDIDATE_MIN_LENGTH,
        max_length=CANDIDATE_MAX_LENGTH,
    )

    # A resumed run rebuilds every derived set, so none is left half-written.
    for name in _candidate_set_dirs(out):
        shutil.rmtree(os.path.join(out, name))
    numbers: Dict[str, int] = {}
    records = []
    for cand in pool.values():
        if cand.kind == _consolidation.ORIGINAL:
            rel = files[cand.parents[0]]
        else:
            set_dir = CANDIDATE_SET_PREFIX + cand.kind
            numbers[cand.kind] = numbers.get(cand.kind, 0) + 1
            rel = os.path.join(set_dir, bin_filename(numbers[cand.kind]))
            ensure_dir(os.path.join(out, set_dir))
            with open(os.path.join(out, rel), "w") as fh:
                fh.writelines(">%s\n%s\n" % (c, sequences[c]) for c in sorted(cand.contigs))
        records.append(
            {
                "file": rel,
                "signature": cand.signature,
                "kind": cand.kind,
                "parents": list(cand.parents),
                "contigs": sorted(cand.contigs),
            }
        )
    work = os.path.join(out, "work_files")
    ensure_dir(work)
    path = os.path.join(work, POOL_FILE)
    with open(path + ".tmp", "w") as fh:
        json.dump({"lengths": dict(sorted(lengths.items())), "candidates": records}, fh)
    os.replace(path + ".tmp", path)
    kinds = _consolidation.count_kinds(pool.values())
    comm(
        "candidate pool: %d candidates from %d input bins (%s)"
        % (len(pool), len(bins), ", ".join("%d %s" % (n, k) for k, n in kinds.items()))
    )


def _read_stats_rows(stats_path: str) -> Dict[str, List[str]]:
    """``{bin name: row}`` from a ``.stats.tsv``, or an empty dict if there is none."""
    rows: Dict[str, List[str]] = {}
    if not os.path.isfile(stats_path):
        return rows
    with open(stats_path) as fh:
        fh.readline()
        for line in fh:
            cut = line.rstrip("\n").split("\t")
            if len(cut) >= 3:
                rows[cut[0]] = cut
    return rows


def _consolidate_candidates(out: str, comp: float, cont: float, weight: float) -> None:
    """Select the default consolidator's bins into ``binsO`` from the scored candidate pool.

    Each candidate's score is read from the ``.stats.tsv`` of the set its file is in. The selected
    bins share no contig, so they are not dereplicated. ``work_files/consolidated_bins.tsv`` records the
    candidate behind each bin, how it was made and the score it was chosen on.
    """
    work = os.path.join(out, "work_files")
    with open(os.path.join(work, POOL_FILE)) as fh:
        saved = json.load(fh)
    lengths = saved["lengths"]
    pool: Dict[str, _consolidation.Candidate] = {}
    files: Dict[str, str] = {}
    stats_rows: Dict[str, Dict[str, List[str]]] = {}
    scores: Dict[str, Tuple[float, float]] = {}
    rows: Dict[str, List[str]] = {}
    for record in saved["candidates"]:
        key = record["signature"]
        pool[key] = _consolidation.Candidate(
            key, frozenset(record["contigs"]), record["kind"], tuple(record["parents"])
        )
        files[key] = record["file"]
        set_dir = os.path.dirname(record["file"])
        if set_dir not in stats_rows:
            stats_rows[set_dir] = _read_stats_rows(os.path.join(out, set_dir + STATS_SUFFIX))
        row = stats_rows[set_dir].get(os.path.splitext(os.path.basename(record["file"]))[0])
        if row is not None:
            scores[key] = (float(row[1]), float(row[2]))
            rows[key] = row
    unscored = len(pool) - len(scores)
    if unscored:
        warning("%d candidate bin(s) have no CheckM score and are not considered" % unscored)

    result = _consolidation.select(
        pool, scores, comp, cont, weight=weight, lookahead=True, lengths=lengths
    )
    binsO = os.path.join(out, "binsO")
    ensure_dir(binsO)
    stats_lines = ["\t".join(_checkm.STATS_HEADER) + "\n"]
    table = ["bin\tkind\tparents\tcompleteness\tcontamination\tscore\tcandidate\n"]
    for number, picked in enumerate(result.chosen, start=1):
        key = picked.candidate.signature
        shutil.copy(os.path.join(out, files[key]), os.path.join(binsO, bin_filename(number)))
        stem = os.path.splitext(bin_filename(number))[0]
        stats_lines.append("\t".join([stem] + rows[key][1:]) + "\n")
        table.append(
            "%s\t%s\t%s\t%s\t%s\t%.4f\t%s\n"
            % (
                stem,
                picked.candidate.kind,
                ",".join(picked.candidate.parents),
                rows[key][1],
                rows[key][2],
                picked.score,
                files[key],
            )
        )
    with open(binsO + STATS_SUFFIX, "w") as fh:
        fh.writelines(stats_lines)
    with open(os.path.join(work, CONSOLIDATED_BINS_FILE), "w") as fh:
        fh.writelines(table)
    kinds = _consolidation.count_kinds(result.chosen)
    comm(
        "selected %d bins (%s) at contamination weight %g; %d candidate(s) overlapped a "
        "better one, %d derived candidate(s) were refused by the lookahead check"
        % (
            len(result.chosen),
            ", ".join("%d %s" % (n, k) for k, n in kinds.items()),
            weight,
            result.blocked,
            result.rejected_by_lookahead,
        )
    )


def _refine(out: str, env: Optional[str], inputs: Sequence[str], dest: str) -> None:
    """Run Binning_refiner on any number of bin folders (relative to *out*) into *dest*.

    Two or three folders are passed as -1/-2/-3, as in earlier versions. More than three use the
    repeatable -i flag.
    """
    refined = "Refined_" + dest[len("bins") :]  # binsAB -> Refined_AB
    if len(inputs) == 2:
        cmd = fill_command(BINNING_REFINER, b1=inputs[0], b2=inputs[1], out=refined)
    elif len(inputs) == 3:
        cmd = fill_command(BINNING_REFINER_3, b1=inputs[0], b2=inputs[1], b3=inputs[2], out=refined)
    else:
        cmd = fill_command(
            BINNING_REFINER_N,
            bins=Raw(" ".join("-i %s" % shlex.quote(b) for b in inputs)),
            out=refined,
        )
    run(cmd, env=None, tool="binning_refiner", cwd=out)
    produced = os.path.join(out, refined, "Refined")
    if dry_run():
        return
    if not os.path.isdir(produced):
        error("Binning_refiner did not produce refined bins for %s." % dest)
    # Replaces any partial set left by a crashed run when --resume redoes the refine step.
    fs_replace_dir(produced, os.path.join(out, dest))
    comm("there are %d refined bins in %s" % (_count_bins(os.path.join(out, dest)), dest))
    shutil.rmtree(os.path.join(out, refined))


def _normalize_bin_extensions(bin_set_dir: str) -> None:
    """Give every bin in a bin-set directory the canonical BIN_EXTENSION.

    One directory can hold ``.fa``, ``.fna`` and ``.fasta``. CheckM takes a single ``-x``
    extension and would score only matching files.
    """
    for f in listdir(bin_set_dir):
        stem, ext = os.path.splitext(f)
        if ext and ext in FASTA_EXTENSIONS and ext != BIN_EXTENSION:
            os.replace(
                os.path.join(bin_set_dir, f), os.path.join(bin_set_dir, stem + BIN_EXTENSION)
            )


def good_bins(stats_path: str, comp: float, cont: float) -> List[str]:
    """Return names of bins in a .stats file meeting the completion/contamination thresholds.

    Both bounds are inclusive. Also used by `metawrap2 score`, so both commands count good bins
    the same way.
    """
    good: List[str] = []
    if not os.path.isfile(stats_path):
        return good
    with open(stats_path) as fh:
        for line in fh:
            if "completeness" in line or "compl" in line:
                continue
            cut = line.rstrip("\n").split("\t")
            if len(cut) < 3:
                continue
            if comp <= float(cut[1]) <= 100 and 0 <= float(cut[2]) <= cont:
                good.append(cut[0])
    return good


#: Memory in GB budgeted per parallel CheckM1 call on the reduced tree (the default with
#: --checkm1) and on the full tree (--full-tree), from CheckM's guidance. pplacer's memory is per
#: process, ~34 GB for the full tree, and it gains little past 8 threads, so sets run in parallel
#: instead of giving one call more threads.
CHECKM1_CALL_GB = _checkm.REDUCED_TREE_GB
CHECKM1_FULL_TREE_CALL_GB = _checkm.FULL_TREE_GB

#: Memory in GB budgeted per parallel CheckM2 call (one DIAMOND block and the prediction model).
CHECKM2_CALL_GB = 24


def checkm_budget(args) -> int:
    """How many CheckM calls fit at once: memory per call first, then one thread per call."""
    if args.checkm2:
        per_call = CHECKM2_CALL_GB
    else:
        per_call = CHECKM1_FULL_TREE_CALL_GB if args.full_tree else CHECKM1_CALL_GB
    return max(1, min(args.memory // per_call, args.threads))


#: The score cache under work_files, one per scorer and CheckM1 tree, so a resume never mixes
#: their scores.
SCORE_CACHE = "checkm_scores.%s.json"

#: Records which scorer filled the per-set ``.stats`` files, so ``--resume`` can refuse to
#: continue a run with another one.
SCORER_FILE = "scorer.txt"

#: Scorer names (see :func:`scorer_name`), as named in messages and as the options that select them.
_SCORER_LABELS = {
    "checkm2": ("CheckM2", "--checkm2 (the default)"),
    "checkm1_reduced": ("CheckM1 on the reduced tree", "--checkm1"),
    "checkm1": ("CheckM1 on the full tree", "--checkm1 --full-tree"),
}


def scorer_name(args) -> str:
    """``checkm2``, ``checkm1_reduced`` or ``checkm1`` (full tree): the scorer *args* select.

    Used in file names. The reduced and full trees can choose different marker sets, so they are
    kept apart. Plain ``checkm1`` is the full tree, which runs from before 2.1.4 used by default.
    """
    if args.checkm2:
        return "checkm2"
    return "checkm1" if args.full_tree else "checkm1_reduced"


def _scorer_of_stats(out: str) -> Optional[str]:
    """The scorer that wrote the per-set ``.stats.tsv`` files in *out*, or None if none exist.

    For runs from before :data:`SCORER_FILE` existed, whose CheckM1 default was the full tree.
    """
    for parent in (out, os.path.join(out, "work_files")):
        for name in sorted(listdir(parent)):
            if not (name.startswith("bins") and name.endswith(STATS_SUFFIX)):
                continue
            found = _checkm.stats_scorer(os.path.join(parent, name))
            if found:
                return found
    return None


def _check_resumed_scorer(out: str, args, resume: bool) -> None:
    """Stop a ``--resume`` that would mix scores from two scorers (or CheckM1 trees) in one run.

    Before 2.1.4 the default was CheckM1, so resuming an older run without ``--checkm1`` would
    re-score only the final set with CheckM2. A run that is not resumed records its scorer.
    """
    path = os.path.join(out, "work_files", SCORER_FILE)
    wanted = scorer_name(args)
    if resume:
        recorded = None
        if os.path.isfile(path):
            with open(path) as fh:
                recorded = fh.read().strip() or None
        else:
            recorded = _scorer_of_stats(out)
        if recorded and recorded != wanted:
            label, flag = _SCORER_LABELS.get(recorded, (recorded, "the same scorer"))
            error(
                "%s was scored with %s, and this run would score with %s.\n"
                "Resume with %s, or start again with --force."
                % (out, label, _SCORER_LABELS[wanted][0], flag)
            )
        if os.path.isfile(path):
            return
    if not dry_run():
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(wanted + "\n")


#: Written under work_files by a default-consolidator run, so ``--resume`` can tell its output from
#: metaWRAP 1's method, whose runs (including all runs from before 2.1.6) never write it.
CONSOLIDATOR_FILE = "consolidator.txt"


def _check_resumed_consolidator(out: str, args, ckpt) -> None:
    """Stop a ``--resume`` that would continue a run with the other consolidation method.

    The two methods build different candidate sets, so continuing one with the other mixes them.
    A metaWRAP 1 run is recognised by its finished Binning_refiner step.
    """
    path = os.path.join(out, "work_files", CONSOLIDATOR_FILE)
    if ckpt.resume:
        recorded = None
        if os.path.isfile(path):
            with open(path) as fh:
                recorded = fh.read().strip() or None
        elif ckpt.finished("refine"):
            recorded = MW1_CONSOLIDATOR
        if recorded and recorded != args.consolidator and args.consolidator != "dastool":
            if recorded == MW1_CONSOLIDATOR:
                error(
                    "%s was started with metaWRAP 1's refinement (the default up to 2.1.5).\n"
                    "Resume with --mw1-refinement, or start again with --force to use the "
                    "current default." % out
                )
            error(
                "%s was started with the default consolidator (metawrap2).\n"
                "Resume without --mw1-refinement, or start again with --force." % out
            )
    if dry_run():
        return
    if args.consolidator != DEFAULT_CONSOLIDATOR:
        if not ckpt.resume and os.path.isfile(path):
            os.remove(path)  # left by a default-consolidator run this one replaces (--force)
        return
    if not ckpt.resume or not os.path.isfile(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(args.consolidator + "\n")


def checkm_jobs(args, n_sets: int) -> int:
    """How many of *n_sets* bin sets to score at once: ``--checkm-jobs``, or the memory budget."""
    asked = getattr(args, "checkm_jobs", 0) or 0
    if asked:
        return max(1, min(asked, n_sets))
    return max(1, min(n_sets, checkm_budget(args)))


def _run_checkm_set(
    bin_set_dir: str,
    args,
    env: Optional[str],
    checkm2_env: Optional[str] = None,
    cache: Optional[_checkm.BinScoreCache] = None,
    threads: Optional[int] = None,
    mem_gb: Optional[int] = None,
) -> str:
    """Run CheckM on one bin set and report the number of 'good' bins. Returns the .stats path.

    Uses `run_checkm_cached` because the hybrid sets overlap: on a three-binner refinement, 32% of
    the bins across the seven sets are repeats. Reported numbers are unchanged.
    """
    comm(
        "Running %s on %s bins"
        % ("CheckM2" if args.checkm2 else "CheckM", os.path.basename(bin_set_dir))
    )
    stats = _checkm.run_checkm_cached(
        bin_set_dir,
        cache=cache,
        prefilter=not args.score_everything,
        say=comm,
        threads=args.threads if threads is None else threads,
        mem_gb=args.memory if mem_gb is None else mem_gb,
        env=env,
        reduced_tree=not args.full_tree,
        checkm2=args.checkm2,
        checkm2_env=checkm2_env,
        ttable=args.checkm2_ttable,
    )
    n = len(good_bins(stats, args.completeness, args.contamination))
    comm(
        "%d 'good' bins found in %s (>%s%% completion and <%s%% contamination)"
        % (n, os.path.basename(bin_set_dir), args.completeness, args.contamination)
    )
    return stats


def _score_bin_sets(
    out: str,
    sets: Sequence[str],
    args,
    env: Optional[str],
    checkm2_env: Optional[str],
    score_cache: _checkm.BinScoreCache,
    jobs: int,
    per_threads: int,
    per_memory: int,
) -> None:
    """Score every bin set, *jobs* at a time. A failure in any set stops the run.

    With *jobs* 1 the sets are scored one after another. Otherwise `plan_scoring` first assigns
    each bin to one set, in set order, so the result does not depend on which set finishes first.
    The sets are then scored in parallel into the cache, and a final pass writes each set's
    ``.stats.tsv`` from the cache without running CheckM.
    """
    if jobs <= 1:
        for name in bar(sets, desc="CheckM on bin sets", unit=" set", total=len(sets)):
            _run_checkm_set(os.path.join(out, name), args, env, checkm2_env, cache=score_cache)
        return

    dirs = [os.path.join(out, name) for name in sets]
    plan = _checkm.plan_scoring(dirs, score_cache, prefilter=not args.score_everything)
    comm(
        "%d bin(s) to score across %d set(s): %s"
        % (
            sum(len(v) for v in plan.values()),
            len(dirs),
            ", ".join("%s %d" % (os.path.basename(d), len(plan[d])) for d in dirs),
        )
    )

    scorer_kwargs = {
        "threads": per_threads,
        "mem_gb": per_memory,
        "env": env,
        "reduced_tree": not args.full_tree,
        "checkm2": args.checkm2,
        "checkm2_env": checkm2_env,
        "ttable": args.checkm2_ttable,
    }

    def score(directory: str) -> None:
        _checkm.score_bins(directory, plan[directory], score_cache, **scorer_kwargs)

    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(score, d): d for d in dirs if plan[d]}
        for future in concurrent.futures.as_completed(futures):
            future.result()  # re-raises a failure from that set
            done += 1
            comm(
                "scored %d of %d bin set(s) (%s)"
                % (done, len(futures), os.path.basename(futures[future]))
            )

    # Every bin is in the cache now, so this writes the per-set files and runs no CheckM.
    for name in sets:
        _run_checkm_set(
            os.path.join(out, name),
            args,
            env,
            checkm2_env,
            cache=score_cache,
            threads=per_threads,
            mem_gb=per_memory,
        )


def _bin_set_dirs(out: str) -> List[str]:
    """Bin-set directories in *out* (binsA, ..., binsAB, ...), sorted, excluding binsM/binsO."""
    dirs = []
    for name in sorted(listdir(out)):
        if (
            name.startswith("bins")
            and name not in ("binsM", "binsO")
            and os.path.isdir(os.path.join(out, name))
        ):
            dirs.append(name)
    return dirs


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    resolve_consolidator(args)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    resolve_memory_gb(args, settings)
    absolutize_paths(args)
    env = env_for("bin_refinement", settings)
    # CheckM2 cannot share an env with CheckM1, so it has its own, looked up separately. The
    # scorer's database is checked here too, before hours of work.
    checkm2_env = scorer_preflight(args, settings) if args.run_checkm else None
    if args.run_checkm:
        scorer_notes(args)
    if args.run_checkm and not args.checkm2 and args.checkm2_ttable:
        comm("--checkm2-ttable applies to CheckM2 only; CheckM1 runs without it.")
    # GUNC and DAS Tool are in this module's env, so the default GUNC step needs no extra install.
    gunc_env = env
    # An external consolidator needs the assembly; check before copying bins or running prodigal.
    # DAS Tool is installed with CheckM1, so its env needs no separate check.
    external = args.consolidator == "dastool"
    candidates = args.consolidator == "metawrap2"
    external_env = env if external else None
    if external and not args.assembly:
        error(
            "--consolidator %s reads the assembly contigs, so it needs -a/--assembly.\n"
            "Pass the assembly the bins were made from. The default consolidator "
            "(--consolidator metawrap2) uses only the bins and does not need -a."
            % args.consolidator
        )
    if external and not args.consolidate:
        # Reported as a conflict instead of resolved by precedence.
        error(
            "--consolidator %s and --skip-consolidation cannot be used together: one picks a "
            "consolidator, the other turns consolidation off.\n"
            "Drop --skip-consolidation to run %s, or drop --consolidator to pick the best single "
            "bin set without mixing bins." % (args.consolidator, args.consolidator)
        )
    given = collect_bin_sets(args)
    if args.consolidator == MW1_CONSOLIDATOR and len(given) > MW1_MAX_BIN_SETS:
        error(
            "%d bin sets given, but --mw1-refinement takes at most %d.\n"
            "metaWRAP 1's method refines and scores combinations of the bin sets with CheckM, so "
            "its cost grows quickly with more sets. The default consolidator (without "
            "--mw1-refinement) takes any number of bin sets. Or refine in groups: consolidate "
            "the first few, then consolidate that result with the rest."
            % (len(given), MW1_MAX_BIN_SETS)
        )
    bin_sets = drop_empty_bin_sets(given)
    if not bin_sets:
        rec = start_run("bin_refinement", args, env, settings, inputs=list(given))
        try:
            if not dry_run():
                final = write_empty_result(args.output, args.completeness, args.contamination)
                comm("0 bins -> %s" % final)
            announcement("BIN_REFINEMENT FINISHED WITH NO BINS")
        finally:
            finish_run(rec, args.output)
        return 0
    validate_inputs(
        "bin_refinement",
        [
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS),
                path,
                "bin set %s" % bin_set_label(i),
            )
            for i, path in enumerate(bin_sets)
        ]
        + ([(check_fasta, args.assembly, "assembly (-a)")] if external else []),
    )
    rec = start_run("bin_refinement", args, env, settings, inputs=list(bin_sets))

    ckpt = make_checkpoint(args.output)
    try:
        # Checked up front because GUNC runs last, after hours of CheckM.
        gunc_db = settings.db("GUNC_DB") if args.run_gunc else ""
        if args.run_gunc and not (gunc_db and os.path.exists(gunc_db)):
            # A dry run warns and keeps printing the plan; a real run stops.
            (warning if dry_run() else error)(
                "GUNC runs by default, but its database is not configured (GUNC_DB=%s). %s"
                % (gunc_db or "<not set>", _gunc.DB_HINT)
            )

        comp, cont = args.completeness, args.contamination
        folders = [(bin_set_label(i), src) for i, src in enumerate(bin_sets)]
        for _label, src in folders:
            if not os.path.isdir(src):
                error("%s is not a valid directory. Exiting." % src)

        announcement("BEGIN PIPELINE!")
        out = args.output
        comm("setting up output folder and copying over bins...")
        # Only under --force. --resume reuses the existing output to skip completed steps; cleaning
        # unconditionally deleted out/bins* on every resume between "refine" and "finalize".
        if os.path.isdir(out) and force():
            warning("%s already exists; removing the old bin sets from it." % out)
            # Remove every bins* working directory and sidecar file (including binsM and binsO), in
            # the output root and in work_files, where a finished run leaves them. Otherwise a rerun
            # with fewer inputs keeps a stale binsC and its binsC.stats.tsv, which the stats tables
            # and plots read.
            for parent in (out, os.path.join(out, "work_files")):
                for name in listdir(parent):
                    path = os.path.join(parent, name)
                    if not name.startswith("bins"):
                        continue
                    if os.path.isdir(path):
                        shutil.rmtree(path, ignore_errors=True)
                    else:
                        os.remove(path)
            # Remove the previous final set and its sidecars, matched by pattern because the name
            # includes the -c/-x thresholds. A leftover set is ambiguous, and copytree in the layout
            # step fails if the directory exists.
            for name in listdir(out):
                if is_refined_bins_dirname(name) and os.path.isdir(os.path.join(out, name)):
                    shutil.rmtree(os.path.join(out, name), ignore_errors=True)
                    for suffix in (STATS_SUFFIX, CONTIGS_SUFFIX):
                        sidecar = os.path.join(out, name + suffix)
                        if os.path.isfile(sidecar):
                            os.remove(sidecar)
            # Remove the copies of input sets under their input folder names. They are only written
            # if absent, so a rerun with changed inputs kept the old ones (a 188-bin MaxBin2 set was
            # shown beside results from the corrected 129-bin set).
            for _label, src in folders:
                base = os.path.basename(strip_trailing_slash(src))
                copy = os.path.join(out, base)
                # Skip if this is the input itself (-A can point inside the output directory).
                if os.path.realpath(copy) == os.path.realpath(src):
                    continue
                if os.path.isdir(copy):
                    shutil.rmtree(copy, ignore_errors=True)
                for suffix in (STATS_SUFFIX, CONTIGS_SUFFIX):
                    sidecar = copy + suffix
                    if os.path.isfile(sidecar):
                        os.remove(sidecar)
        ensure_dir(out)

        n_binnings = 0
        for dest, src in folders:
            if not os.path.isdir(src):
                error("%s is not a valid directory. Exiting." % src)
            kept = _copy_bins(src, os.path.join(out, dest))
            comm("there are %d bins in %s" % (kept, dest))
            if kept == 0 and not dry_run():
                error("No valid bin sets were given.")
            n_binnings += 1
        comm("There are %d bin sets!" % n_binnings)

        # One CheckM score cache for the whole refinement: the hybrid sets overlap (32% identical
        # bins across the seven sets of a three-binner run). Kept under work_files for --resume.
        # One file per scorer, so a resume with the other scorer cannot reuse these scores.
        # Bins scored by earlier runs with the same scorer, version and database are reused from
        # the shared cache (metawrap2.scorecache), unless --no-score-cache.
        shared_scores = (
            open_score_cache(
                args,
                settings,
                scorer_name(args),
                ttable=args.checkm2_ttable if args.checkm2 else None,
                env=(checkm2_env or env) if args.checkm2 else env,
            )
            if args.run_checkm
            else None
        )
        score_cache = _checkm.BinScoreCache(
            os.path.join(out, "work_files", SCORE_CACHE % scorer_name(args)), shared=shared_scores
        )
        if args.run_checkm:
            _check_resumed_scorer(out, args, ckpt.resume)
        _check_resumed_consolidator(out, args, ckpt)

        comm("Removing special characters from contig names...")
        rewritten = 0
        for dest, _src in folders:
            for f in listdir(os.path.join(out, dest)):
                if _fix_contig_naming(os.path.join(out, dest, f)):
                    rewritten += 1
        comm("normalised contig headers in %d bin file(s); the rest were already clean" % rewritten)

        if external:
            comm(
                "Skipping bin refinement: %s consolidates the input bin sets directly, without "
                "building combinations." % args.consolidator
            )
        elif candidates and n_binnings >= 2:
            if dry_run():
                comm("(dry run) would build the candidate pool from the %d bin sets" % n_binnings)
            elif ckpt.todo("candidate_pool"):
                announcement("BUILDING CANDIDATE BINS")
                _build_candidate_pool(out, [label for label, _src in folders])
                ckpt.done("candidate_pool")
            else:
                comm("skipping the candidate pool (already built; --resume)")
        elif args.refine and n_binnings >= 2:
            if ckpt.todo("refine"):
                announcement("BEGIN BIN REFINEMENT")
                plan = refinement_plan(
                    [label for label, _src in folders], combinations=args.refine_combinations
                )
                comm(
                    "%d bin folders give %d combination(s) to refine (%s). "
                    "Refining each of them."
                    % (
                        n_binnings,
                        len(plan),
                        ", ".join(
                            "+".join(lbl[len("bins") :] for lbl in inputs) for inputs, _dest in plan
                        ),
                    )
                )
                for inputs, dest in bar(
                    plan, desc="refining bin-set combinations", unit=" combo", total=len(plan)
                ):
                    _refine(out, env, inputs, dest)
                comm("Bin refinement finished successfully!")
                ckpt.done("refine")
            else:
                comm("skipping bin refinement (already done; --resume)")
        elif not args.refine:
            comm("Skipping bin refinement; continuing with the %d bin sets given." % n_binnings)
        else:
            comm("There is only one bin folder, so there is nothing to refine.")

        comm("renaming bin files to the .fa extension...")
        for name in _bin_set_dirs(out):
            _normalize_bin_extensions(os.path.join(out, name))

        comm("checking that every refined bin set contains bins...")
        for name in [] if dry_run() else _bin_set_dirs(out):
            if _count_bins(os.path.join(out, name)) == 0:
                comm("Removing bin set %s because it has 0 refined bins" % name)
                shutil.rmtree(os.path.join(out, name))

        # ── CheckM on all bin sets ──────────────────────────────────────────────────────────
        if external:
            # Not run: CheckM on every input set only serves the built-in consolidation, and an
            # external consolidator scores bins itself. It would add the module's largest cost.
            comm(
                "Skipping CheckM on the input bin sets: %s scores the bins itself. The "
                "consolidated set is scored with CheckM below, as with the default consolidator."
                % args.consolidator
            )
        elif args.run_checkm:
            if ckpt.todo("checkm"):
                announcement("RUNNING CHECKM ON ALL SETS OF BINS")
                # CheckM on each of up to 7 independent bin sets, several at once when memory allows.
                sets = _bin_set_dirs(out) or ([d for d, _ in folders] if dry_run() else [])
                jobs = checkm_jobs(args, len(sets)) if sets else 1
                per_threads = max(1, args.threads // jobs)
                per_memory = max(1, args.memory // jobs)
                if jobs > 1:
                    comm(
                        "Scoring %d bin sets %d at a time, with %d thread(s) and %d GB each "
                        "(--checkm-jobs 1 scores them one at a time)"
                        % (len(sets), jobs, per_threads, per_memory)
                    )
                _score_bin_sets(
                    out,
                    sets,
                    args,
                    env,
                    checkm2_env,
                    score_cache,
                    jobs,
                    per_threads,
                    per_memory,
                )
                ckpt.done("checkm")
            else:
                comm("skipping CheckM on bin sets (already done; --resume)")
        else:
            comm("Skipping CheckM. Bins cannot be consolidated without CheckM scores.")

        # ── consolidate ─────────────────────────────────────────────────────────────────────
        best_bin_set: Optional[str] = None
        if external:
            if ckpt.todo("consolidate"):
                announcement("CONSOLIDATING ALL BIN SETS WITH " + args.consolidator.upper())
                _run_dastool(args, external_env, out, [label for label, _src in folders])
                ckpt.done("consolidate")
            else:
                comm("skipping %s consolidation (already done; --resume)" % args.consolidator)
            best_bin_set = "binsO"
        elif candidates and n_binnings >= 2:
            announcement(
                "CONSOLIDATING ALL BIN SETS FROM INTERSECTION, UNION AND DIFFERENCE CANDIDATES"
            )
            if dry_run():
                comm(
                    "(dry run) would select non-overlapping candidates into binsO by completeness "
                    "- %g * contamination" % args.contamination_weight
                )
            elif ckpt.todo("consolidate"):
                binsO = os.path.join(out, "binsO")
                # A crashed attempt leaves a partial binsO.
                shutil.rmtree(binsO, ignore_errors=True)
                _consolidate_candidates(out, comp, cont, args.contamination_weight)
                ckpt.done("consolidate")
            else:
                comm("skipping consolidation (already done; --resume)")
            best_bin_set = "binsO"
        elif args.consolidate and args.run_checkm:
            announcement("CONSOLIDATING ALL BIN SETS BY CHOOSING THE BEST VERSION OF EACH BIN")
            if n_binnings == 1:
                comm("There is only one original bin folder, so there is nothing to consolidate.")
                best_bin_set = "binsA"
            else:
                if ckpt.todo("consolidate") and not dry_run():
                    comm("Consolidating %d original bin folders and the refined bins." % n_binnings)
                    binsM = os.path.join(out, "binsM")
                    # Clear binsM (and binsO below) left by a crashed consolidation on --resume.
                    shutil.rmtree(binsM, ignore_errors=True)
                    shutil.copytree(os.path.join(out, "binsA"), binsM)
                    shutil.copy(os.path.join(out, "binsA" + STATS_SUFFIX), binsM + STATS_SUFFIX)
                    for name in _bin_set_dirs(out):
                        if name == "binsA":
                            continue
                        stats = os.path.join(out, name + STATS_SUFFIX)
                        if not os.path.isfile(stats):
                            continue
                        comm("merging %s and binsM" % name)
                        tmp = os.path.join(out, "binsM1")
                        _refinement.consolidate(
                            binsM,
                            os.path.join(out, name),
                            binsM + STATS_SUFFIX,
                            stats,
                            tmp,
                            comp,
                            cont,
                        )
                        shutil.rmtree(binsM)
                        os.remove(binsM + STATS_SUFFIX)
                        os.replace(tmp, binsM)
                        os.replace(tmp + STATS_SUFFIX, binsM + STATS_SUFFIX)

                    binsO = os.path.join(out, "binsO")
                    # A crashed attempt leaves a partial binsO; dereplicate writes into it, so stale
                    # bins would reach the final set, and the "keep" move would fail. Cleared once
                    # for all three modes.
                    shutil.rmtree(binsO, ignore_errors=True)
                    if args.dereplicate == "keep":
                        comm("Skipping dereplication of contigs between bins (--keep-ambiguous)")
                        os.replace(binsM, binsO)
                        os.replace(binsM + STATS_SUFFIX, binsO + STATS_SUFFIX)
                    elif args.dereplicate == "partial":
                        comm(
                            "Finding contigs in more than one bin and keeping each only in "
                            "the best bin..."
                        )
                        _refinement.dereplicate(binsM + STATS_SUFFIX, binsM, binsO, mode="best")
                        shutil.copy(binsM + STATS_SUFFIX, binsO + STATS_SUFFIX)
                    elif args.dereplicate == "remove":
                        comm(
                            "Finding contigs in more than one bin and removing them from all "
                            "bins..."
                        )
                        _refinement.dereplicate(binsM + STATS_SUFFIX, binsM, binsO, mode="remove")
                        shutil.copy(binsM + STATS_SUFFIX, binsO + STATS_SUFFIX)
                    ckpt.done("consolidate")
                elif dry_run():
                    comm(
                        "(dry run) would consolidate the bin sets into binsM, then "
                        "dereplicate contigs into binsO"
                    )
                else:
                    comm("skipping bin consolidation (already done; --resume)")
                best_bin_set = "binsO"
        elif not args.consolidate:
            comm("Skipping bin consolidation; picking the best bin set without mixing bins.")
            if not args.run_checkm:
                comm("CheckM was not run, so the best bin set cannot be chosen; using binsA.")
                best_bin_set = "binsA"
            else:
                best, max_good = None, -1
                for name in _bin_set_dirs(out):
                    n = len(good_bins(os.path.join(out, name + STATS_SUFFIX), comp, cont))
                    comm("%d 'good' bins found in %s" % (n, name))
                    if n > max_good:
                        max_good, best = n, name
                best_bin_set = best
                comm("the best bin set (most good bins) is %s" % best_bin_set)
        else:
            best_bin_set = "binsA"
        comm("The best non-reassembled versions of the bins are in %s" % best_bin_set)

        # ── finalize: re-QC the consolidated set and lay out the final outputs ───────────────
        if ckpt.todo("finalize") and not dry_run():
            # With the built-in consolidation, --keep-ambiguous skips the re-score because binsO is
            # binsM unchanged and already scored. An external consolidator's output has not been
            # scored, so it is always re-scored. The default's bins are scored candidates, so its
            # re-score comes from the score cache and runs no CheckM.
            rescore = args.run_checkm and best_bin_set == "binsO"
            if rescore and args.consolidator == MW1_CONSOLIDATOR and args.dereplicate == "keep":
                rescore = False
            if rescore:
                binsO = os.path.join(out, "binsO")
                if _checkm.count_bins(binsO) == 0:
                    error(
                        "Consolidation produced no bins: no bin from any bin set met "
                        "the quality thresholds of >%s%% completeness and <%s%% "
                        "contamination.\n"
                        "The per-bin-set scores are in %s/bins*.stats. Check the best "
                        "completeness there, and re-run with a lower -c (or a higher "
                        "-x) if that suits your analysis. If even the best completeness "
                        "is very low, the assembly or coverage is usually too thin to "
                        "recover genomes." % (comp, cont, out)
                    )
                comm("Re-running CheckM on binsO bins")
                stats = _checkm.run_checkm_cached(
                    binsO,
                    cache=score_cache,
                    prefilter=not args.score_everything,
                    say=comm,
                    threads=args.threads,
                    mem_gb=args.memory,
                    env=env,
                    reduced_tree=not args.full_tree,
                    checkm2=args.checkm2,
                    checkm2_env=checkm2_env,
                    ttable=args.checkm2_ttable,
                )
                comm("Removing bins below the quality thresholds...")
                good = set(good_bins(stats, comp, cont))
                for f in list(os.listdir(binsO)):
                    if f.endswith(BIN_EXTENSION) and os.path.splitext(f)[0] not in good:
                        print(
                            "Removing %s: it is below the quality thresholds after "
                            "contig dereplication" % os.path.splitext(f)[0]
                        )
                        os.remove(os.path.join(binsO, f))
                # rewrite stats to keep only the good bins
                with open(stats) as fh:
                    lines = fh.readlines()
                with open(stats, "w") as fh:
                    fh.write(lines[0])
                    for line in lines[1:]:
                        cut = line.rstrip("\n").split("\t")
                        if (
                            len(cut) >= 3
                            and comp <= float(cut[1]) <= 100
                            and 0 <= float(cut[2]) <= cont
                        ):
                            fh.write(line)
                comm(
                    "Bin quality re-scored after contig dereplication: %d bins pass the "
                    "thresholds." % len(good)
                )

            if args.run_checkm:
                comm("making completion and contamination ranking plots of final outputs")
                # Runs with cwd=out (the plot is written to the working directory), so paths must
                # be absolute; a relative -o resolved to out/out/binsA.stats.
                stats_files = quoted_words(
                    sorted(
                        os.path.abspath(os.path.join(out, f))
                        for f in os.listdir(out)
                        if f.endswith(STATS_SUFFIX)
                    )
                )
                ensure_dir(os.path.join(out, "figures"))
                run(
                    fill_command(PLOT_BINNING, comp=int(comp), cont=int(cont), stats=stats_files),
                    env=None,
                    tool="plot_binning_results",
                    cwd=out,
                )
                png = os.path.join(out, "binning_results.png")
                if os.path.isfile(png):
                    os.replace(png, os.path.join(out, "figures", "binning_results.png"))

            # ── move intermediate files aside and lay out the final outputs ──────────────────
            announcement("MOVING OVER TEMPORARY FILES")
            if n_binnings != 1:
                work = os.path.join(out, "work_files")
                ensure_dir(work)
                # Every working bin set and sidecar, by prefix (a fixed name list missed binsD/binsE
                # with more than three inputs).
                for f in list(os.listdir(out)):
                    if f.startswith("bins"):
                        # The destination exists if this folder was refined before; a plain move
                        # raised FileExistsError at the end of the run.
                        fs_replace_dir(os.path.join(out, f), os.path.join(work, f))

            final = os.path.join(out, refined_bins_dirname(comp, cont))
            if best_bin_set is None:
                # Every branch above sets this; checked for the type checker and future branches.
                error("internal error: no best bin set was chosen")
            best_dir = (
                os.path.join(out, "work_files", best_bin_set)
                if n_binnings != 1
                else os.path.join(out, best_bin_set)
            )
            if os.path.isdir(best_dir):
                # With two or more sets, consolidation already applied the thresholds. With one set
                # they are applied here; previously this folder received every input bin.
                single_set_stats = best_dir + STATS_SUFFIX
                if n_binnings == 1 and args.run_checkm and os.path.isfile(single_set_stats):
                    kept = _refinement.select_passing(best_dir, single_set_stats, final, comp, cont)
                    comm(
                        "%d of the %d bins pass %d%% completeness and %d%% contamination"
                        % (kept, len(fasta_files(best_dir)), int(comp), int(cont))
                    )
                else:
                    shutil.copytree(best_dir, final)
                    if os.path.isfile(single_set_stats):
                        shutil.copy(single_set_stats, final + STATS_SUFFIX)

            # restore the original input bin sets under their input folder names
            for dest, src in folders:
                base = os.path.basename(strip_trailing_slash(src))
                src_dir = (
                    os.path.join(out, "work_files", dest)
                    if n_binnings != 1
                    else os.path.join(out, dest)
                )
                if os.path.isdir(src_dir) and not os.path.exists(os.path.join(out, base)):
                    shutil.copytree(src_dir, os.path.join(out, base))
                    if os.path.isfile(src_dir + STATS_SUFFIX):
                        shutil.copy(src_dir + STATS_SUFFIX, os.path.join(out, base + STATS_SUFFIX))

            if args.run_checkm:
                comm("making contig membership files (for Anvio and other applications)")
                for name in os.listdir(out):
                    d = os.path.join(out, name)
                    if name.endswith("_bins") and os.path.isdir(d):
                        write_contig_table(d, d + CONTIGS_SUFFIX)
            ckpt.done("finalize")
        elif dry_run():
            if args.run_checkm:
                ensure_dir(os.path.join(out, "figures"))
                run(
                    fill_command(
                        PLOT_BINNING,
                        comp=int(comp),
                        cont=int(cont),
                        stats=os.path.join(out, "<bin-set>.stats"),
                    ),
                    env=None,
                    tool="plot_binning_results",
                    cwd=out,
                )
            comm("(dry run) would lay out the final bins in %s" % refined_bins_dirname(comp, cont))
        else:
            comm("skipping final output layout (already done; --resume)")

        # ── GUNC: chimerism, which CheckM's marker genes do not detect ──────────────────────
        # A separate checkpoint from "finalize", so a resumed run can add GUNC without redoing the
        # layout. Chimeric bins are annotated and reported, not removed.
        final_bins = os.path.join(out, refined_bins_dirname(comp, cont))
        if args.run_gunc and ckpt.todo("gunc"):
            announcement("ASSESSING THE FINAL BINS FOR CHIMERISM WITH GUNC")
            if dry_run():
                # A placeholder so the printed command shows where the database goes.
                _gunc.run_gunc(
                    final_bins,
                    threads=args.threads,
                    db=gunc_db or "<GUNC_DB>",
                    env=gunc_env,
                    out_dir=os.path.join(out, "work_files", "gunc"),
                )
            elif _checkm.count_bins(final_bins) == 0:
                comm("no final bins to assess; skipping GUNC")
                ckpt.done("gunc")
            else:
                report = _gunc.run_gunc(
                    final_bins,
                    threads=args.threads,
                    db=gunc_db,
                    env=gunc_env,
                    out_dir=os.path.join(out, "work_files", "gunc"),
                )
                verdicts = _gunc.parse_report(str(report))
                chimeric = sorted(name for name, v in verdicts.items() if not v.passed)
                n = _gunc.annotate_stats(final_bins + STATS_SUFFIX, verdicts)
                comm(
                    "GUNC scored %d bins; %d annotated in %s"
                    % (len(verdicts), n, os.path.basename(final_bins) + STATS_SUFFIX)
                )
                if chimeric:
                    warning(
                        "GUNC calls %d of %d final bins chimeric (clade separation score above "
                        "%s): %s\nThey are kept, with gunc_pass=no in the stats table. A bin that "
                        "spans two genomes can still score as complete and clean in CheckM."
                        % (len(chimeric), len(verdicts), _gunc.CSS_THRESHOLD, ", ".join(chimeric))
                    )
                else:
                    comm("GUNC found no chimeric bins among the final %d" % len(verdicts))
                ckpt.done("gunc")
        elif args.run_gunc:
            comm("skipping GUNC (already done; --resume)")

        # ── MIMAG tiers for the final bin set ───────────────────────────────────────────────
        # Tagged here, not in checkm.run_checkm, so only the final set gets tiers and not the
        # candidate hybrid sets.
        #
        # Refinement runs before annotate, so no 16S search has been done. Bins passing 90/5 are
        # tagged high_quality_unverified_rrna (completeness passed, rRNA not measured) instead of
        # "medium". Re-running after `annotate --bakta` resolves them.
        if not dry_run() and os.path.isdir(final_bins):
            stats_path = final_bins + STATS_SUFFIX
            assessments = _checkm.tag_mimag(stats_path)
            if assessments:
                comm(_mimag.summarise(assessments))
                counts = _mimag.tally(assessments)
                if counts[_mimag.MimagTag.HIGH_UNVERIFIED]:
                    comm(_mimag.NO_EVIDENCE_HINT)

        # Recovery and chimerism caveats for the final set (see metawrap2.recovery).
        if not dry_run() and args.assembly and os.path.isdir(final_bins):
            recovery = _recovery.measure(args.assembly, final_bins)
            table = os.path.join(args.output, "recovery" + TSV_SUFFIX)
            with open(table, "w", newline="") as fh:
                writer = csv.DictWriter(
                    fh,
                    fieldnames=_recovery.RECOVERY_COLUMNS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow(_recovery.recovery_row(recovery))
            comm(recovery.describe())
            comm(_recovery.UNBINNED_CAVEAT)
            comm("-> %s" % table)
        comm(_recovery.CHIMERISM_CAVEAT)

        announcement("BIN_REFINEMENT PIPELINE FINISHED SUCCESSFULLY!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
