# Verification evidence

Reproduced by `python -m pytest tests -q` (27 tests) and
`python examples/run_verification.py` (figures + GIF + CSVs/JSON under
`docs/results/`). Verified on Windows 11, Python 3.13.5, NumPy 2.1.3;
CI on Ubuntu.

## What each test establishes

| Test group | Property established |
| --- | --- |
| Series-wall hand values | Two-layer wall reproduces the hand-computed R = 0.05 K/W, q = 1600 W, and interface temperature 68 degC exactly (1e-12). |
| Matrix vs direct reduction | The node/matrix solution equals the independent series/parallel reduction (temperatures AND heat rate) on the recorded wall. |
| Single segment | No-interior-solve path works (q = 4000 W hand case). |
| Energy balance | Wall dissipates exactly what enters: hot-side supply = cold-side removal = q; every interior node closes (residual ~1e-12 W measured on the recorded run). |
| Parallel superposition | Parallel branches carry hand-computed shares (200/400 W of 600 W) that sum exactly to the wall rate; conductances superpose as 1/(1/3)+1/(1/6) = 6. |
| Parallel helper | `parallel_equivalent_resistance([0.5, 0.25])` matches the product formula. |
| Parallel wall balance | Parallel-group wall closes energy balance to 1e-12. |
| Fin analytic core | (in the fin tests within the same suite) midpoint symmetry θ(L/2) from base = from tip exactly (0.0 K measured); L→0 limit recovers plain base convection q = h·Ac·θ0; η ∈ (0,1]; ε > 1 for the recorded case; monotone temperature; insulated/corrected/convective heat rates ordered consistently. |
| Validation | NaN/inf/non-numeric/zero/negative physical inputs rejected with specific messages; parallel groups need ≥ 2 branches; nested parallel rejected; positions outside [0, L] rejected; hot ≤ cold rejected. |

## Recorded run (`docs/results/`)

- `wall_profile.png/.svg` + `wall_nodes.csv` + `wall_links.csv`:
  five-element wall (films + brick + contact + insulation) at 100→20 degC:
  q = 70.97 W, energy-balance residual **1.7e-12 W**.
- `fin_solution.png/.svg`: aluminum fin (50×3 mm, L = 150 mm, k = 200,
  h = 60): insulated/corrected/convective-tip temperature profiles and the
  spanwise heat decay; efficiency 0.4427; midpoint symmetry error 0.0 K.
- `fin_sankey.png` + `fin_sankey.csv`: balanced 4-segment energy Sankey —
  inflow 38.328436 W = Σ outflows, residual **0.0 W** (computed model
  output, one consistent unit).
- `fin_length_progression.gif`: 41 frames, each a full analytic solution
  at that fin length (parameter progression, not decorative).
- `run_summary.json`: the metrics above + runtime context.

## Limits of this evidence

- Analytic/1-D: no comparison with multi-dimensional or measured data.
- The Sankey shows computed model output, not measured energy.
- Fin limits: Bi ≪ 1 assumed; tip-correction is an engineering
  approximation whose accuracy vs the exact convective tip is itself
  quantified in the recorded run.
