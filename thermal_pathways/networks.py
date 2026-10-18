"""Series/parallel thermal resistance networks with node/matrix assembly.

Formulation (1-D, steady state, no internal heat generation)::

    conduction through a slab :  R = L / (k A)     [K/W]
    convection at a surface   :  R = 1 / (h A)     [K/W]
    contact resistance        :  R_c given directly [K/W]
    series chain              :  R_tot = sum(R_i)
    parallel paths            :  1/R_eq = sum(1/R_j)   (conductances add)

Node method: every network node ``k`` obeys ``sum_j G_kj (T_k - T_j) = 0``
with conductance ``G = 1/R``.  The two boundary node temperatures are known
(the hot-side and cold-side reservoirs); the interior node temperatures are
obtained from a small dense linear solve (:func:`numpy.linalg.solve`) after
assembling the conductance matrix.  For the networks constructible here
(series chains plus at most one parallel group between shared flanking
nodes) the same answer follows from direct series/parallel reduction; the
test suite cross-checks the matrix solution against that reduction and
verifies the energy balance at every node to roundoff.

Parallel branches are each reduced to their own series resistance and then
connected between the same two flanking nodes; the heat rates of the
branches therefore add exactly to the heat rate through the group.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple

import numpy as np

from .validation import require_finite, require_positive

__all__ = [
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
]


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ConductionSegment:
    """One-dimensional conduction through a slab, ``R = L / (k A)`` [K/W].

    Parameters
    ----------
    name : str
        Label used in results and figures.
    thickness : float
        Layer thickness ``L`` [m], strictly positive.
    conductivity : float
        Thermal conductivity ``k`` [W/(m K)], strictly positive.
    area : float
        Cross-sectional (heat-flow) area ``A`` [m2], strictly positive.
    """

    name: str
    thickness: float
    conductivity: float
    area: float

    def __post_init__(self) -> None:
        require_positive(f"{self.name}.thickness", self.thickness)
        require_positive(f"{self.name}.conductivity", self.conductivity)
        require_positive(f"{self.name}.area", self.area)
        object.__setattr__(self, "thickness", float(self.thickness))
        object.__setattr__(self, "conductivity", float(self.conductivity))
        object.__setattr__(self, "area", float(self.area))

    @property
    def resistance(self) -> float:
        """Thermal resistance ``L / (k A)`` [K/W]."""
        return self.thickness / (self.conductivity * self.area)


@dataclass(frozen=True)
class ConvectionSegment:
    """Convection at a surface, ``R = 1 / (h A)`` [K/W]."""

    name: str
    h: float
    area: float

    def __post_init__(self) -> None:
        require_positive(f"{self.name}.h", self.h)
        require_positive(f"{self.name}.area", self.area)
        object.__setattr__(self, "h", float(self.h))
        object.__setattr__(self, "area", float(self.area))

    @property
    def resistance(self) -> float:
        """Film resistance ``1 / (h A)`` [K/W]."""
        return 1.0 / (self.h * self.area)


@dataclass(frozen=True)
class ContactSegment:
    """An interface (contact) resistance specified directly in [K/W].

    ``resistance`` is the *total* resistance of the joint for the whole area
    it applies to (not per unit area).
    """

    name: str
    resistance: float

    def __post_init__(self) -> None:
        require_positive(f"{self.name}.resistance", self.resistance)
        object.__setattr__(self, "resistance", float(self.resistance))


@dataclass(frozen=True)
class SeriesBranch:
    """One branch of a parallel group: a series chain of plain segments.

    The branch is reduced to a single resistance ``sum(R_i)``; nested
    parallel groups are not allowed (documented 1-D restriction).
    """

    name: str
    segments: Tuple[object, ...]

    def __post_init__(self) -> None:
        if not self.segments:
            raise ValueError(f"branch '{self.name}' must contain at least one segment")
        for segment in self.segments:
            if isinstance(segment, ParallelSegment):
                raise ValueError(
                    f"branch '{self.name}': nested ParallelSegment is not supported"
                )
            if not hasattr(segment, "resistance"):
                raise TypeError(
                    f"branch '{self.name}': segment {segment!r} has no resistance"
                )
        object.__setattr__(self, "segments", tuple(self.segments))

    @property
    def resistance(self) -> float:
        """Series resistance of the branch [K/W]."""
        return float(sum(segment.resistance for segment in self.segments))


@dataclass(frozen=True)
class ParallelSegment:
    """Parallel heat-flow paths sharing two flanking nodes.

    ``1/R_eq = sum(1/R_branch)``.  Example: insulation and timber studs
    occupying different areas of the same wall cavity.
    """

    name: str
    branches: Tuple[SeriesBranch, ...]

    def __post_init__(self) -> None:
        if len(self.branches) < 2:
            raise ValueError(
                f"parallel group '{self.name}' needs at least two branches"
            )
        for branch in self.branches:
            if not isinstance(branch, SeriesBranch):
                raise TypeError(
                    f"parallel group '{self.name}': branches must be SeriesBranch"
                )
        object.__setattr__(self, "branches", tuple(self.branches))

    @property
    def resistance(self) -> float:
        """Equivalent resistance of the parallel group [K/W]."""
        return parallel_equivalent_resistance(
            [branch.resistance for branch in self.branches]
        )


def series_resistance(segments: Sequence[object]) -> float:
    """Total resistance of a series chain: ``sum(R_i)`` [K/W]."""
    return float(sum(segment.resistance for segment in segments))


def parallel_equivalent_resistance(resistances: Sequence[float]) -> float:
    """Equivalent resistance of parallel paths: ``1 / sum(1/R_j)`` [K/W]."""
    values = [require_positive("resistance", value) for value in resistances]
    conductance = sum(1.0 / value for value in values)
    return 1.0 / conductance


# ---------------------------------------------------------------------------
# Solution containers
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LinkResult:
    """A single two-terminal element (or parallel branch) between two nodes."""

    name: str
    resistance: float          # K/W
    node_a: int                # hotter-side node index in chain order
    node_b: int                # colder-side node index in chain order
    heat_rate: float           # W, positive from node_a to node_b
    temperature_drop: float    # K, T(node_a) - T(node_b)


@dataclass(frozen=True)
class WallSolution:
    """Steady-state solution of a composite-wall resistance network.

    Attributes
    ----------
    node_labels : tuple of str
        Labels of all nodes in chain order (boundaries first and last).
    node_temperatures : numpy.ndarray
        Temperatures [degC] of all nodes, same order as ``node_labels``.
    links : tuple of LinkResult
        Every element/branch with its resistance, heat rate and drop.
    total_resistance : float
        Series/parallel-reduced resistance of the whole wall [K/W].
    heat_rate : float
        Steady heat rate through the wall [W].
    """

    node_labels: Tuple[str, ...]
    node_temperatures: np.ndarray
    links: Tuple[LinkResult, ...]
    total_resistance: float
    heat_rate: float
    hot_side_temperature: float
    cold_side_temperature: float
    _adjacency: dict = field(default_factory=dict, repr=False, compare=False)

    def heat_rates_for_group(self, group_prefix: str) -> Tuple[LinkResult, ...]:
        """All link results whose name starts with ``group_prefix``."""
        return tuple(link for link in self.links if link.name.startswith(group_prefix))

    def hot_side_supply(self) -> float:
        """Heat rate leaving the hot boundary node into the wall [W]."""
        return float(
            sum(
                conductance * (self.node_temperatures[0] - self.node_temperatures[j])
                for j, conductance in self._adjacency.get(0, ())
            )
        )

    def cold_side_removal(self) -> float:
        """Heat rate arriving at the cold boundary node from the wall [W]."""
        last = len(self.node_temperatures) - 1
        return float(
            sum(
                conductance * (self.node_temperatures[j] - self.node_temperatures[last])
                for j, conductance in self._adjacency.get(last, ())
            )
        )

    def max_interior_node_imbalance(self) -> float:
        """Max |sum of fluxes| over interior nodes [W]; roundoff if converged."""
        worst = 0.0
        for k in range(1, len(self.node_temperatures) - 1):
            imbalance = sum(
                conductance * (self.node_temperatures[k] - self.node_temperatures[j])
                for j, conductance in self._adjacency.get(k, ())
            )
            worst = max(worst, abs(imbalance))
        return float(worst)

    def energy_balance_residual(self) -> float:
        """Wall energy-balance residual [W].

        ``max(|hot-side supply - q|, |cold-side removal - q|, max interior
        node imbalance)``: the network must dissipate exactly what enters,
        to roundoff.
        """
        return float(
            max(
                abs(self.hot_side_supply() - self.heat_rate),
                abs(self.cold_side_removal() - self.heat_rate),
                self.max_interior_node_imbalance(),
            )
        )


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
def solve_wall(
    segments: Sequence[object],
    hot_side_temperature: float,
    cold_side_temperature: float,
    *,
    hot_label: str = "hot-side reservoir",
    cold_label: str = "cold-side reservoir",
) -> WallSolution:
    """Solve a composite-wall network built from ``segments`` in chain order.

    The segments are traversed from the hot side to the cold side.  A
    :class:`ParallelSegment` contributes one link per branch between the
    same two flanking nodes; every other segment contributes one link.
    Interior node temperatures come from the node/matrix assembly described
    in the module docstring (dense linear solve); the total resistance and
    wall heat rate come from the independent series/parallel reduction.

    Parameters
    ----------
    segments : sequence
        Objects exposing ``.resistance`` (see module classes), at least one.
    hot_side_temperature, cold_side_temperature : float
        Boundary temperatures [degC] (any finite values; hot > cold expected).
    hot_label, cold_label : str, optional
        Labels for the boundary nodes.

    Returns
    -------
    WallSolution
    """
    segments = tuple(segments)
    if not segments:
        raise ValueError("segments must contain at least one element")
    t_hot = require_finite("hot_side_temperature", hot_side_temperature)
    t_cold = require_finite("cold_side_temperature", cold_side_temperature)
    if t_hot <= t_cold:
        raise ValueError(
            f"hot_side_temperature ({t_hot}) must exceed cold_side_temperature ({t_cold})"
        )
    for segment in segments:
        if not hasattr(segment, "resistance"):
            raise TypeError(f"segment {segment!r} has no resistance")

    n_nodes = len(segments) + 1
    node_labels = [hot_label]
    for i in range(1, n_nodes - 1):
        node_labels.append(f"interface {i}: {segments[i - 1].name} / {segments[i].name}")
    node_labels.append(cold_label)

    # Expand segments into two-terminal links (conductances) between nodes.
    links: list = []
    for i, segment in enumerate(segments):
        if isinstance(segment, ParallelSegment):
            for branch in segment.branches:
                links.append(
                    (f"{segment.name}: {branch.name}", branch.resistance, i, i + 1)
                )
        else:
            links.append((segment.name, float(segment.resistance), i, i + 1))

    # Total resistance via independent series/parallel reduction.
    total = series_resistance(segments)

    temperatures = np.empty(n_nodes, dtype=float)
    temperatures[0] = t_hot
    temperatures[-1] = t_cold
    interior = n_nodes - 2  # number of unknown node temperatures
    if interior > 0:
        matrix = np.zeros((interior, interior), dtype=float)
        rhs = np.zeros(interior, dtype=float)
        for name, resistance, a, b in links:
            conductance = 1.0 / resistance
            for (k, j) in ((a, b), (b, a)):
                if 1 <= k <= interior:
                    matrix[k - 1, k - 1] += conductance
                    if 1 <= j <= interior:
                        matrix[k - 1, j - 1] -= conductance
                    else:
                        rhs[k - 1] += conductance * temperatures[j]
        temperatures[1:-1] = np.linalg.solve(matrix, rhs)

    heat_rate = (t_hot - t_cold) / total

    adjacency: dict = {}
    for _, resistance, a, b in links:
        conductance = 1.0 / resistance
        adjacency.setdefault(a, []).append((b, conductance))
        adjacency.setdefault(b, []).append((a, conductance))

    link_results = tuple(
        LinkResult(
            name=name,
            resistance=resistance,
            node_a=a,
            node_b=b,
            heat_rate=(temperatures[a] - temperatures[b]) / resistance,
            temperature_drop=temperatures[a] - temperatures[b],
        )
        for name, resistance, a, b in links
    )

    return WallSolution(
        node_labels=tuple(node_labels),
        node_temperatures=temperatures,
        links=link_results,
        total_resistance=total,
        heat_rate=heat_rate,
        hot_side_temperature=t_hot,
        cold_side_temperature=t_cold,
        _adjacency=adjacency,
    )
