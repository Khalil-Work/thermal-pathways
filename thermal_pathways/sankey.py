"""Quantitatively balanced energy Sankey for the straight-fin case.

The fin input heat rate ``q_f`` at the base is split into the convective
losses of consecutive spanwise segments plus the heat reaching the tip.
With the analytic insulated-tip solution the split is exact:

    q(x) = sqrt(h P k Ac) theta0 sinh(m (L - x)) / cosh(m L)

    loss over [x_i, x_i+1] = q(x_i) - q(x_i+1)
    tip flow              = q(L) = 0            (insulated tip)

so ``sum(losses) + tip == q_f`` to floating-point roundoff.  The residual
is computed and recorded on the figure and in ``run_summary.json``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

import numpy as np

from .fins import StraightRectangularFin

__all__ = ["SankeyFlow", "fin_sankey", "draw_sankey"]

PALETTE = ["#2b6a99", "#c8553d", "#4d8b4f", "#8a6d3b", "#6b5b8a", "#5b5b5b"]


@dataclass(frozen=True)
class SankeyFlow:
    """A balanced set of energy flows in one consistent unit (watts)."""

    source_label: str
    inflow: float
    outflows: Tuple[Tuple[str, float], ...]
    residual: float
    unit: str = "W"

    @property
    def outflow_total(self) -> float:
        return float(sum(watts for _, watts in self.outflows))


def fin_sankey(fin: StraightRectangularFin, n_segments: int = 4) -> SankeyFlow:
    """Build the balanced flow set for ``fin`` (insulated tip).

    Raises
    ------
    AssertionError
        If the balance residual exceeds 1e-9 * q_f (it is roundoff by
        construction, so this guards against regressions).
    """
    if n_segments < 1:
        raise ValueError("n_segments must be at least 1")
    edges = np.linspace(0.0, fin.length, n_segments + 1)
    q_at = np.atleast_1d(fin.section_heat_rate(edges))
    outflows: List[Tuple[str, float]] = []
    for lo, hi, qa, qb in zip(edges[:-1], edges[1:], q_at[:-1], q_at[1:]):
        label = f"convection {lo * 1e3:.0f}\u2013{hi * 1e3:.0f} mm"
        outflows.append((label, float(qa - qb)))
    outflows.append(("insulated tip (q(L) = 0)", float(q_at[-1])))

    inflow = fin.heat_rate()
    residual = inflow - sum(watts for _, watts in outflows)
    if abs(residual) > 1e-9 * max(abs(inflow), 1.0):
        raise AssertionError(
            f"Sankey flows do not balance: residual = {residual:.3e} W"
        )
    return SankeyFlow(
        source_label="fin base input q_f",
        inflow=inflow,
        outflows=tuple(outflows),
        residual=residual,
    )


def draw_sankey(flow: SankeyFlow, path: str, *, title: str, dpi: int = 180) -> None:
    """Draw ``flow`` as a left-to-right Sankey diagram (PNG, white background).

    Ribbon and node heights are proportional to the flow in watts — the
    diagram is a quantitative representation of the computed model output,
    not a sketch.
    """
    import matplotlib.pyplot as plt

    watts_per_unit = flow.inflow / 0.78  # axes-fraction height of full inflow
    gap = 0.030
    n = len(flow.outflows)

    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.axis("off")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)

    # Right nodes: stacked top to bottom.
    heights = [max(watts / watts_per_unit, 0.0) for _, watts in flow.outflows]
    total_stack = sum(heights) + gap * (n - 1)
    top = 0.89
    right_bars = []  # (y_bottom, y_top, watts)
    y = top
    for h in heights:
        right_bars.append((y - h, y, h))
        y -= h + gap
    stack_bottom = y + gap

    # Left node spans the whole inflow.
    left_top = top
    left_bottom = top - flow.inflow / watts_per_unit

    x_left_bar = (0.06, 0.115)
    x_right_bar = (0.70, 0.755)

    # Ribbons: subdivide the left bar in the same top-down order.
    y_cursor = left_top
    for (label, watts), (rb, rt, h), color in zip(flow.outflows, right_bars, PALETTE):
        band_height = watts / watts_per_unit
        la2, la1 = y_cursor, y_cursor - band_height
        y_cursor = la1
        ribbon = plt.Polygon(
            [
                (x_left_bar[1], la1),
                (x_left_bar[1], la2),
                (x_right_bar[0], rt),
                (x_right_bar[0], rb),
            ],
            closed=True,
            facecolor=color,
            alpha=0.45,
            edgecolor="none",
        )
        ax.add_patch(ribbon)
        ax.add_patch(
            plt.Rectangle(
                (x_right_bar[0], rb),
                x_right_bar[1] - x_right_bar[0],
                max(rt - rb, 0.002),
                facecolor=color,
                edgecolor="white",
                linewidth=0.5,
            )
        )
        value_note = f"{watts:.3f} W" if watts > 1e-9 else "0.000 W"
        ax.text(
            x_right_bar[1] + 0.012,
            (rb + rt) / 2,
            f"{label}\n{value_note}",
            va="center",
            ha="left",
            fontsize=9,
            color="#222222",
        )

    # Left node bar.
    ax.add_patch(
        plt.Rectangle(
            (x_left_bar[0], left_bottom),
            x_left_bar[1] - x_left_bar[0],
            left_top - left_bottom,
            facecolor="#2b6a99",
            edgecolor="white",
            linewidth=0.5,
        )
    )
    ax.text(
        x_left_bar[0] - 0.012,
        (left_top + left_bottom) / 2,
        f"{flow.source_label}\n{flow.inflow:.3f} W",
        ha="right",
        va="center",
        fontsize=10,
        fontweight="bold",
        color="#222222",
    )

    ax.text(
        0.5,
        0.965,
        title,
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
    )
    ax.text(
        0.5,
        min(stack_bottom, left_bottom) - 0.055,
        (
            f"Node totals balance to roundoff: in {flow.inflow:.6f} W = "
            f"sum(out) {flow.outflow_total:.6f} W  "
            f"(residual {flow.residual:+.2e} W) \u2014 computed model output"
        ),
        ha="center",
        va="center",
        fontsize=9,
        color="#333333",
    )

    fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
    plt.close(fig)
