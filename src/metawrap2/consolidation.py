"""Bin consolidation from union, difference and intersection candidates, selected greedily.

The algorithm behind bin_refinement's default consolidator (``--consolidator metawrap2``).
metaWRAP 1's method (``--mw1-refinement``) builds intersections of the input bin sets with
Binning_refiner, so every candidate is a subset of an
input bin and no output bin can be more complete than the best input version of it. This one adds
unions, which can, and differences:

1. Bins are nodes of an overlap graph. Two bins are adjacent when they share at least
   :data:`MIN_OVERLAP_FRACTION` of the smaller bin's contigs and their lengths are within
   :data:`MAX_LENGTH_RATIO` of each other.
2. Every clique of 2 to :data:`MAX_COMBO` bins gives an intersection, a union and one difference
   per member. Cliques, because a connected component can be a chain (A overlaps B, B overlaps C,
   A and C share nothing), and the union of A and C is two genomes.
3. Derived candidates outside :data:`MIN_LENGTH`..:data:`MAX_LENGTH` bp are dropped; the input bins
   are kept whatever their length. Identical contig sets are kept once, under the first kind that
   produced them (input bins first).
4. Every candidate is scored with CheckM, and candidates are taken best first by
   ``completeness - weight * contamination``, skipping any that shares a contig with one already
   taken. Ties are broken by N50 (higher first), then input bins before derived ones, then the
   candidate's signature.
5. Lookahead guard: before a derived candidate is taken, the candidates it would block are packed
   greedily over input bins and intersections only, and the derived candidate is taken only if it
   scores at least that packing. Without it, the union of A and B in the chain above blocks A and
   C together, scoring less than the two of them.

The candidate pool borrows from two other consolidators. Intersection, union and difference
candidates over cliques of the overlap graph, the length window, pruning of identical contig sets
and the N50 tiebreak follow Binette (Mainguy and Hoede, 2024). Merging overlapping bins into new
candidates also follows MAGScoT (Ruehlemann et al., 2022). The lookahead guard is MetaWrap2's own.

Every collection that is iterated is sorted first, so the result does not depend on
``PYTHONHASHSEED`` or on the order bins were read in.
"""

from __future__ import annotations

import hashlib
import itertools
from dataclasses import dataclass
from typing import Dict, FrozenSet, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

from .checkm import MAX_PLAUSIBLE_BIN_BP, MIN_PLAUSIBLE_BIN_BP
from .constants import CONTAMINATION_WEIGHT
from .refinement import passes

# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: Set operations that build derived candidates.
OPERATIONS = ("intersection", "union", "difference")
#: Two bins are adjacent in the overlap graph when they share at least this fraction of the smaller
#: bin's contigs. With no floor (any shared contig, as in Binette) the pool, and so the CheckM cost,
#: grows several-fold on large samples, and the floor also gave more clean high-quality MAGs.
MIN_OVERLAP_FRACTION = 0.10
#: Two bins whose lengths differ more than this many times are not adjacent, so a large bin is not
#: merged with a small fragment. 0 turns the check off.
MAX_LENGTH_RATIO = 4.0
#: Largest clique combined. Three binners give at most three-way cliques.
MAX_COMBO = 3
#: Length window (bp) for derived candidates, the bounds outside which a bin is not scored.
MIN_LENGTH = MIN_PLAUSIBLE_BIN_BP
MAX_LENGTH = MAX_PLAUSIBLE_BIN_BP
#: Characters of the SHA-256 hex digest kept as a candidate's signature.
SIGNATURE_LENGTH = 24
# ───────────────────────────────────────────────────────────────────────────────────────────

ORIGINAL = "original"
#: Kinds the lookahead guard packs its alternative over: what metaWRAP 1's method can reach.
BASELINE_KINDS = (ORIGINAL, "intersection")


@dataclass(frozen=True)
class Candidate:
    """One candidate bin: its contigs, how it was made, and from which input bins."""

    signature: str
    contigs: FrozenSet[str]
    kind: str
    parents: Tuple[str, ...]


@dataclass(frozen=True)
class Selected:
    """A candidate chosen for the consolidated set, with the scores it was chosen on."""

    candidate: Candidate
    completeness: float
    contamination: float
    score: float


@dataclass
class Selection:
    """The chosen candidates, best first, and counts of what was passed over."""

    chosen: List[Selected]
    blocked: int = 0
    rejected_by_lookahead: int = 0


def signature(contigs) -> str:
    """A stable identifier for a contig set, so identical candidates are kept once."""
    digest = hashlib.sha256("\n".join(sorted(contigs)).encode()).hexdigest()
    return digest[:SIGNATURE_LENGTH]


def total_length(contigs, lengths: Mapping[str, int]) -> int:
    return sum(lengths.get(c, 0) for c in contigs)


def n50(contigs, lengths: Mapping[str, int]) -> int:
    """N50 of a contig set, the first tiebreak between equally scored candidates."""
    sizes = sorted((lengths.get(c, 0) for c in contigs), reverse=True)
    total = sum(sizes)
    if total == 0:
        return 0
    half, running = total / 2.0, 0
    for size in sizes:
        running += size
        if running >= half:
            return size
    return sizes[-1]


def overlap_graph(
    bins: Mapping[str, FrozenSet[str]],
    lengths: Mapping[str, int],
    min_overlap: float = MIN_OVERLAP_FRACTION,
    max_length_ratio: float = MAX_LENGTH_RATIO,
) -> Dict[str, Set[str]]:
    """Adjacency of the bins that share enough contigs to be combined.

    Built from a contig-to-bins index, so pairs that share nothing are never compared.
    """
    by_contig: Dict[str, List[str]] = {}
    for name in sorted(bins):
        for contig in bins[name]:
            by_contig.setdefault(contig, []).append(name)

    shared: Dict[Tuple[str, str], int] = {}
    for names in by_contig.values():
        if len(names) < 2:
            continue
        for pair in itertools.combinations(sorted(set(names)), 2):
            shared[pair] = shared.get(pair, 0) + 1

    adjacency: Dict[str, Set[str]] = {}
    for (a, b), n_shared in sorted(shared.items()):
        smaller = min(len(bins[a]), len(bins[b]))
        if min_overlap > 0 and (not smaller or n_shared / smaller < min_overlap):
            continue
        if max_length_ratio > 0:
            len_a, len_b = total_length(bins[a], lengths), total_length(bins[b], lengths)
            low, high = min(len_a, len_b), max(len_a, len_b)
            if low <= 0 or high / low > max_length_ratio:
                continue
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    return adjacency


def components(adjacency: Mapping[str, Set[str]], names: Sequence[str]) -> List[List[str]]:
    """Connected components of the overlap graph, each sorted, ordered by their first member."""
    seen: Set[str] = set()
    found: List[List[str]] = []
    for start in sorted(names):
        if start in seen:
            continue
        stack, members = [start], []
        seen.add(start)
        while stack:
            node = stack.pop()
            members.append(node)
            for neighbour in sorted(adjacency.get(node, ())):
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        found.append(sorted(members))
    return found


def cliques(
    members: Sequence[str], adjacency: Mapping[str, Set[str]], size: int
) -> Iterator[Tuple[str, ...]]:
    """Every clique of *size* among sorted *members*, in lexicographic order.

    The same tuples as filtering ``itertools.combinations(members, size)`` for cliques, without
    enumerating the combinations that are not.
    """

    def extend(clique: Tuple[str, ...], pool: Sequence[str]) -> Iterator[Tuple[str, ...]]:
        if len(clique) == size:
            yield clique
            return
        for i, node in enumerate(pool):
            neighbours = adjacency.get(node, set())
            yield from extend(clique + (node,), [m for m in pool[i + 1 :] if m in neighbours])

    yield from extend((), sorted(members))


def build_candidates(
    bins: Mapping[str, FrozenSet[str]],
    lengths: Mapping[str, int],
    operations: Sequence[str] = OPERATIONS,
    max_combo: int = MAX_COMBO,
    min_overlap: float = MIN_OVERLAP_FRACTION,
    max_length_ratio: float = MAX_LENGTH_RATIO,
    min_length: int = MIN_LENGTH,
    max_length: int = MAX_LENGTH,
) -> Dict[str, Candidate]:
    """The candidate pool, ``{signature: Candidate}`` in a fixed order: input bins first.

    *bins* maps a bin name to its contig ids, *lengths* a contig id to its length.
    """
    unknown = sorted(set(operations) - set(OPERATIONS))
    if unknown:
        raise ValueError("unknown operation(s): %s" % ", ".join(unknown))
    pool: Dict[str, Candidate] = {}

    def add(contigs: FrozenSet[str], kind: str, parents: Tuple[str, ...]) -> None:
        if not contigs:
            return
        if kind != ORIGINAL:
            length = total_length(contigs, lengths)
            if length < min_length or length > max_length:
                return
        key = signature(contigs)
        if key not in pool:
            pool[key] = Candidate(key, contigs, kind, parents)

    for name in sorted(bins):
        add(frozenset(bins[name]), ORIGINAL, (name,))

    adjacency = overlap_graph(bins, lengths, min_overlap, max_length_ratio)
    for members in components(adjacency, list(bins)):
        for size in range(2, min(max_combo, len(members)) + 1):
            for combo in cliques(members, adjacency, size):
                sets = [frozenset(bins[n]) for n in combo]
                if "intersection" in operations:
                    add(frozenset.intersection(*sets), "intersection", combo)
                if "union" in operations:
                    add(frozenset.union(*sets), "union", combo)
                if "difference" in operations:
                    for i, name in enumerate(combo):
                        others = [s for j, s in enumerate(sets) if j != i]
                        add(sets[i].difference(*others), "difference", (name,) + combo)
    return pool


def select(
    pool: Mapping[str, Candidate],
    scores: Mapping[str, Tuple[float, float]],
    min_completeness: float,
    max_contamination: float,
    weight: float = CONTAMINATION_WEIGHT,
    lookahead: bool = True,
    lengths: Optional[Mapping[str, int]] = None,
) -> Selection:
    """Choose disjoint candidates best first. *scores* maps a signature to (completeness,
    contamination); candidates with no score or failing the thresholds are not considered.
    """
    ranked = []
    for key in sorted(pool):
        if key not in scores:
            continue
        completeness, contamination = scores[key]
        if not passes(completeness, contamination, min_completeness, max_contamination):
            continue
        cand = pool[key]
        score = completeness - weight * contamination
        ranked.append((score, completeness, contamination, cand))
    ranked.sort(
        key=lambda item: (
            -item[0],
            -n50(item[3].contigs, lengths) if lengths else 0,
            0 if item[3].kind == ORIGINAL else 1,
            item[3].signature,
        )
    )

    holders: Dict[str, List[int]] = {}
    for index, item in enumerate(ranked):
        for contig in item[3].contigs:
            holders.setdefault(contig, []).append(index)

    def alternative(index: int, contigs: FrozenSet[str], used: Set[str]) -> float:
        """What the baseline kinds could still score from the contigs this candidate takes.

        Packed over input bins and intersections only, the pool metaWRAP 1's method searches.
        Packing over every kind lets a competing union stand in for the alternative, and then the
        guard never rejects anything.
        """
        blocked: Set[int] = set()
        for contig in contigs:
            blocked.update(holders.get(contig, ()))
        blocked.discard(index)
        total = 0.0
        taken: Set[str] = set()
        for j in sorted(blocked):
            other = ranked[j][3]
            if other.kind not in BASELINE_KINDS:
                continue
            if other.contigs & used or other.contigs & taken:
                continue
            total += ranked[j][0]
            taken |= other.contigs
        return total

    result = Selection(chosen=[])
    used: Set[str] = set()
    for index, (score, completeness, contamination, cand) in enumerate(ranked):
        if cand.contigs & used:
            result.blocked += 1
            continue
        if lookahead and cand.kind != ORIGINAL and score < alternative(index, cand.contigs, used):
            result.rejected_by_lookahead += 1
            continue
        result.chosen.append(Selected(cand, completeness, contamination, score))
        used |= cand.contigs
    return result


def count_kinds(candidates) -> Dict[str, int]:
    """``{kind: count}`` over Candidates or Selected, in OPERATIONS order after the input bins."""
    counts = {kind: 0 for kind in (ORIGINAL,) + OPERATIONS}
    for item in candidates:
        cand = item.candidate if isinstance(item, Selected) else item
        counts[cand.kind] = counts.get(cand.kind, 0) + 1
    return counts
