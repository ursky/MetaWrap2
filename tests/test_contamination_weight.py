"""One contamination weight for every bin ranking.

``choose_best_bin`` (reassemble_bins) and :func:`metawrap2.utils.quality_score` (bin_refinement's
consolidation and dereplication) rank bins by completeness against weighted contamination. Both
read the weight from :mod:`metawrap2.constants`, so changing it there changes both.
"""

from __future__ import annotations

import itertools
import os
import random
import subprocess
import sys
from pathlib import Path

from metawrap2 import constants, utils
from metawrap2.scripts import choose_best_bin as cbb


def test_choose_best_bin_uses_the_constant_from_constants():
    assert cbb.CONTAMINATION_WEIGHT is constants.CONTAMINATION_WEIGHT
    assert utils.CONTAMINATION_WEIGHT is constants.CONTAMINATION_WEIGHT


def test_changing_the_constant_changes_choose_best_bin(monkeypatch):
    clean = cbb.BinVersion("bin_1", "orig", 80.0, 1.0, 1000)
    fuller = cbb.BinVersion("bin_1", "strict", 85.0, 2.5, 1000)
    assert clean.beats(fuller)  # 1.5 more contamination outweighs 5 more completeness at 5x
    monkeypatch.setattr(cbb, "CONTAMINATION_WEIGHT", 2.0)
    assert fuller.beats(clean)


def test_both_scores_order_bins_the_same():
    """completeness + W*(100 - contamination) is quality_score plus 100*W."""
    rng = random.Random(7)
    versions = [
        cbb.BinVersion("b", str(i), round(rng.uniform(50, 100), 2), round(rng.uniform(0, 10), 2), 1)
        for i in range(200)
    ]
    for a, b in itertools.combinations(versions, 2):
        by_script = (a.score > b.score) - (a.score < b.score)
        qa = utils.quality_score(a.completeness, a.contamination)
        qb = utils.quality_score(b.completeness, b.contamination)
        by_utils = (qa > qb) - (qa < qb)
        if abs(qa - qb) > 1e-9:
            assert by_script == by_utils, (a, b)


def test_the_script_still_runs_as_a_module(tmp_path):
    """`python -m metawrap2.scripts.choose_best_bin`, as reassemble_bins calls it."""
    stats = tmp_path / "s.stats.tsv"
    stats.write_text(
        "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
        "bin_1.orig\t80.0\t1.0\t0.5\tx\t1000\t2000\n"
        "bin_1.strict\t85.0\t2.5\t0.5\tx\t1000\t2000\n"
    )
    out = subprocess.run(
        [sys.executable, "-m", "metawrap2.scripts.choose_best_bin", str(stats), "70", "10"],
        capture_output=True,
        text=True,
        check=True,
        env=dict(os.environ, PYTHONPATH=str(Path(cbb.__file__).parents[2])),
    )
    assert out.stdout.split() == ["bin_1.orig"]
