"""Reproducible verification run for Thermal Pathways.

Cases:
1. Composite wall (two layers + contact resistance): temperature profile,
   hand-verifiable heat rate, roundoff energy balance.
2. Straight rectangular fin (documented aluminum case): temperature
   distribution, heat rates (insulated/corrected/convective tip),
   efficiency/effectiveness, midpoint symmetry.
3. Balanced fin Sankey (4 segments): node totals close to roundoff.

Outputs under ``docs/results/``: wall profile figure (cover), fin figure,
Sankey PNG, fin GIF (parameter progression over fin length - each frame a
full analytic solution), CSVs, run_summary.json.
"""

from __future__ import annotations

import json
import pathlib
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as mpl_animation
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from thermal_pathways import (
    ConductionSegment,
    ContactSegment,
    ConvectionSegment,
    ParallelSegment,
    SeriesBranch,
    StraightRectangularFin,
    draw_sankey,
    fin_sankey,
    solve_wall,
)

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "docs" / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

BLUE, RED, GREEN, GRAY = "#2b6a99", "#c8553d", "#4d8b4f", "#555555"
plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
    }
)


def wall_case():
    segments = [
        ConvectionSegment("hot film", h=40.0, area=2.0),
        ConductionSegment("brick", thickness=0.10, conductivity=0.7, area=2.0),
        ContactSegment("mortar joint", resistance=0.01),
        ConductionSegment("insulation", thickness=0.08, conductivity=0.04, area=2.0),
        ConvectionSegment("cold film", h=15.0, area=2.0),
    ]
    return segments, solve_wall(segments, 100.0, 20.0)


def fin_case(length=0.15):
    return StraightRectangularFin(
        length=length, width=0.05, thickness=0.003,
        conductivity=200.0, h=60.0,
        base_temperature=120.0, fluid_temperature=30.0,
    )


def make_wall_figure(segments, solution):
    # Cumulative thickness for the x-axis (contact/convection = zero width).
    xs, labels = [0.0], ["base"]
    for seg in segments:
        width = getattr(seg, "thickness", 0.0)
        xs.append(xs[-1] + width)
        labels.append(seg.name)
    temps = solution.node_temperatures

    fig, ax = plt.subplots(figsize=(9.6, 4.4))
    ax.plot(xs, temps, "o-", color=BLUE, lw=2.0, ms=6)
    for x, T, label in zip(xs, temps, labels):
        ax.annotate(f"{T:.2f} degC", (x, T), xytext=(6, 6),
                    textcoords="offset points", fontsize=8, color="#333333")
    ax.set_xlabel("cumulative wall thickness [m] (films/contact shown at their node)")
    ax.set_ylabel("temperature [degC]")
    ax.set_title(f"Composite-wall temperature profile "
                 f"(q = {solution.heat_rate:.2f} W)")
    fig.tight_layout()
    fig.savefig(RESULTS / "wall_profile.png", dpi=150, bbox_inches="tight")
    fig.savefig(RESULTS / "wall_profile.svg", dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_fin_figure(fin):
    x = np.linspace(0.0, fin.length, 241)
    fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.2))
    ax[0].plot(x * 1e3, fin.temperature(x), "-", color=RED, lw=1.8,
               label="insulated tip (exact)")
    ax[0].plot(x * 1e3, fin.temperature_excess_convective_tip(x) + fin.fluid_temperature,
               "--", color=GRAY, lw=1.4, label="convective tip (exact)")
    ax[0].axhline(fin.fluid_temperature, color=GREEN, ls=":", lw=1.2,
                  label="fluid")
    ax[0].set_xlabel("distance from base [mm]")
    ax[0].set_ylabel("T [degC]")
    ax[0].set_title("Fin temperature distribution")
    ax[0].legend(frameon=False, fontsize=8)

    q = fin.section_heat_rate(x)
    ax[1].plot(x * 1e3, q, "-", color=BLUE, lw=1.8)
    ax[1].set_xlabel("distance from base [mm]")
    ax[1].set_ylabel("conductive heat rate q(x) [W]")
    ax[1].set_title(f"Spanwise heat decay (q(0) = {q[0]:.3f} W)")
    fig.tight_layout()
    fig.savefig(RESULTS / "fin_solution.png", dpi=150, bbox_inches="tight")
    fig.savefig(RESULTS / "fin_solution.svg", dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_gif(frames=41):
    lengths = np.linspace(0.005, 0.30, frames)
    fig, ax = plt.subplots(figsize=(7.6, 4.8), dpi=140)
    writer = mpl_animation.PillowWriter(fps=8)
    with writer.saving(fig, RESULTS / "fin_length_progression.gif", dpi=140):
        for L in lengths:
            fin = fin_case(length=L)
            x = np.linspace(0.0, L, 120)
            ax.clear()
            ax.plot(x * 1e3, fin.temperature(x), "-", color=RED, lw=1.8)
            ax.axhline(fin.fluid_temperature, color=GREEN, ls=":", lw=1.0)
            ax.set_xlim(0.0, 300.0)
            ax.set_ylim(29.0, 121.0)
            ax.set_xlabel("distance from base [mm]")
            ax.set_ylabel("T [degC]")
            ax.set_title(f"Fin length L = {L * 1e3:.0f} mm: "
                         f"q_f = {fin.heat_rate():.3f} W, "
                         f"eta = {fin.efficiency():.3f}")
            writer.grab_frame()
    plt.close(fig)
    return frames


def write_csvs(solution, fin, flow):
    rows = ["node_index,label,temperature_degC"]
    for i, (label, T) in enumerate(zip(solution.node_labels,
                                       solution.node_temperatures)):
        rows.append(f"{i},{label},{T:.9f}")
    (RESULTS / "wall_nodes.csv").write_text("\n".join(rows), encoding="utf-8")

    rows = ["element,resistance_K_per_W,heat_rate_W,temperature_drop_K"]
    for link in solution.links:
        rows.append(f"{link.name},{link.resistance:.9f},{link.heat_rate:.9f},"
                    f"{link.temperature_drop:.9f}")
    (RESULTS / "wall_links.csv").write_text("\n".join(rows), encoding="utf-8")

    rows = ["label,watts"]
    for label, watts in flow.outflows:
        rows.append(f"{label},{watts:.9f}")
    rows.append(f"inflow,{flow.inflow:.9f}")
    rows.append(f"residual,{flow.residual:+.6e}")
    (RESULTS / "fin_sankey.csv").write_text("\n".join(rows), encoding="utf-8")


def main():
    segments, solution = wall_case()
    make_wall_figure(segments, solution)

    fin = fin_case()
    make_fin_figure(fin)
    flow = fin_sankey(fin, n_segments=4)
    draw_sankey(flow, str(RESULTS / "fin_sankey.png"),
                title="Straight-fin energy balance (computed model output)")
    frames = make_gif()
    write_csvs(solution, fin, flow)

    # Independent hand-check of the two-layer core (films/contact excluded).
    core = sum(seg.resistance for seg in segments[1:4])
    summary = {
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "system": platform.platform(),
        },
        "wall": {
            "total_resistance_K_per_W": solution.total_resistance,
            "heat_rate_W": solution.heat_rate,
            "energy_balance_residual_W": solution.energy_balance_residual(),
            "core_series_resistance_K_per_W": core,
        },
        "fin": {
            "case": "aluminum 50x3 mm, L=150 mm, k=200, h=60, Tb=120, Tinf=30",
            "m_per_m": fin.m,
            "q_insulated_W": fin.heat_rate(),
            "q_corrected_W": fin.heat_rate_corrected(),
            "q_convective_tip_W": fin.heat_rate_convective_tip(),
            "efficiency": fin.efficiency(),
            "effectiveness": fin.effectiveness(),
            "midpoint_symmetry_K": abs(
                fin.temperature_excess(fin.length / 2)
                - fin.temperature_excess_from_tip(fin.length / 2)
            ),
        },
        "sankey": {
            "inflow_W": flow.inflow,
            "outflow_total_W": flow.outflow_total,
            "residual_W": flow.residual,
            "note": "computed model output; one consistent quantity [W]",
        },
        "animation": {"file": "fin_length_progression.gif", "frames": frames,
                      "content": "each frame is a full analytic solution at "
                                  "that fin length"},
    }
    (RESULTS / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("wall", "fin", "sankey")}, indent=2))
    print("frames", frames)


if __name__ == "__main__":
    main()
