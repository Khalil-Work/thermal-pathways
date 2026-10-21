# Model and method

## Composite-wall resistance networks (1-D, steady state)

Conduction through a slab: R = L/(kA) [K/W]; convection at a surface:
R = 1/(hA); contact resistance given directly. Series chains add
resistances; parallel paths between shared flanking nodes add
conductances (1/R_eq = Σ 1/R_j). Nested parallel groups are excluded by
design (documented 1-D restriction).

**Node method:** every interior node k obeys Σ_j G_kj (T_k − T_j) = 0.
The two boundary temperatures are known; interior temperatures come from
a small dense linear solve. The total resistance and heat rate come from
the *independent* direct series/parallel reduction — the test suite
cross-checks the two methods against each other and verifies energy
balance at every node to roundoff (measured 1.7e-12 W on the recorded
case).

## Straight rectangular fin (exact analytic)

With perimeter P = 2w + 2t, cross-section Ac = wt, and fin parameter
m = sqrt(hP/(kAc)):

- Insulated tip: θ(x) = θ0 cosh(m(L−x))/cosh(mL);
  q_f = sqrt(h P k Ac) θ0 tanh(mL).
- Corrected length Lc = L + Ac/P (tip convection via the adiabatic
  algebra): q_f,c = sqrt(hPkAc) θ0 tanh(mLc); η_f = tanh(mLc)/(mLc);
  ε_f = q_f,c/(h Ac θ0).
- Exact convective-tip solution (H = h/(mk)):
  θ(x) = θ0 [cosh(m(L−x)) + H sinh(m(L−x))]/[cosh(mL) + H sinh(mL)].

Assumptions: 1-D (Bi = ht/k ≪ 1), steady, constant properties, no
radiation. Section heat rate q(x) = sqrt(hPkAc) θ0 sinh(m(L−x))/cosh(mL)
decreases monotonically to zero at the insulated tip.

## Balanced fin Sankey

The base input q_f splits into consecutive spanwise convective losses
plus the tip flow; with the analytic section rates the split is exact, so
Σ(losses) + tip = q_f to roundoff (residual computed, printed on the
figure, recorded in run_summary.json). Ribbon and node heights are
proportional to watts — one consistent quantity, computed model output
(not measured energy).

## Limits

1-D steady conduction/convection only; no radiation, no transient, no
temperature-dependent properties, no fins other than straight
rectangular. Networks: series + one-level parallel only.
