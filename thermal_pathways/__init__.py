"""thermal_pathways — resistance networks, fins and balanced heat flows.

Original Python/NumPy implementation of

* series/parallel composite-wall thermal resistance networks with
  node/matrix assembly and roundoff energy balances
  (:mod:`thermal_pathways.networks`),
* the straight rectangular fin with exact analytic solutions
  (:mod:`thermal_pathways.fins`),
* a quantitatively balanced energy Sankey for the fin case
  (:mod:`thermal_pathways.sankey`).
"""

from .fins import StraightRectangularFin
from .networks import (
    ConductionSegment,
    ContactSegment,
    ConvectionSegment,
    LinkResult,
    ParallelSegment,
    SeriesBranch,
    WallSolution,
    parallel_equivalent_resistance,
    series_resistance,
    solve_wall,
)
from .sankey import SankeyFlow, draw_sankey, fin_sankey
from .validation import require_finite, require_positive

__version__ = "1.0.0"

__all__ = [
    "StraightRectangularFin",
    "ConductionSegment",
    "ConvectionSegment",
    "ContactSegment",
    "SeriesBranch",
    "ParallelSegment",
    "LinkResult",
    "WallSolution",
    "solve_wall",
    "series_resistance",
    "parallel_equivalent_resistance",
    "SankeyFlow",
    "fin_sankey",
    "draw_sankey",
    "require_positive",
    "require_finite",
    "__version__",
]
