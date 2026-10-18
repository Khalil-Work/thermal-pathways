# API reference

Package `thermal_pathways` (modules re-exported at package level).

## Networks (`thermal_pathways.networks`)

- `ConductionSegment(name, thickness, conductivity, area)` — R = L/(kA) [K/W].
- `ConvectionSegment(name, h, area)` — R = 1/(hA).
- `ContactSegment(name, resistance)` — direct [K/W].
- `SeriesBranch(name, segments)` — series chain (no nested parallel).
- `ParallelSegment(name, branches)` — parallel paths, ≥ 2 SeriesBranch.
- `solve_wall(segments, hot_side_temperature, cold_side_temperature)` →
  `WallSolution` with `node_temperatures`, `links` (per element/branch:
  resistance, heat rate, drop), `total_resistance`, `heat_rate`,
  `hot_side_supply()`, `cold_side_removal()`, `max_interior_node_imbalance()`,
  `energy_balance_residual()` — all in W.
- `series_resistance(segments)`, `parallel_equivalent_resistance(list)`.

## Fins (`thermal_pathways.fins`)

`StraightRectangularFin(length, width, thickness, conductivity, h,
base_temperature, fluid_temperature)` with:
- `temperature(x)` / `temperature_excess(x)` (insulated tip);
- `temperature_excess_convective_tip(x)`; `temperature_excess_from_tip(y)`;
- `heat_rate()` / `heat_rate_corrected()` / `heat_rate_convective_tip()`;
- `section_heat_rate(x)`; `efficiency()`; `effectiveness()`;
- properties `perimeter`, `cross_section_area`, `fin_parameter` (m),
  `corrected_length`.

## Sankey (`thermal_pathways.sankey`)

- `fin_sankey(fin, n_segments=4)` → `SankeyFlow` (raises if the balance
  residual exceeds 1e-9·q_f — it is roundoff by construction).
- `draw_sankey(flow, path, title=...)` — white-background PNG with the
  residual printed on the figure.

## Validation (`thermal_pathways.validation`)

`require_positive`, `require_finite` — NaN/inf/non-numeric rejected
everywhere; temperatures may be negative (degC).

## Example (executed)

```python
from thermal_pathways import ConductionSegment, solve_wall

wall = [ConductionSegment("brick", 0.10, 0.7, 2.0),
        ConductionSegment("insulation", 0.08, 0.04, 2.0)]
s = solve_wall(wall, 100.0, 20.0)
print(s.heat_rate, s.energy_balance_residual())
```
