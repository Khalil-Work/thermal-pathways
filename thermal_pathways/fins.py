"""Straight rectangular fin — exact analytic solutions (no Bessel functions).

Geometry and conventions
------------------------
A straight rectangular fin of length ``L`` [m], width ``w`` [m] and
thickness ``t`` [m] extends from a base at ``T_b`` into a fluid at
``T_inf`` with uniform convection coefficient ``h``::

    perimeter              P  = 2w + 2t          [m]
    cross-section area     Ac = w * t            [m2]
    fin parameter          m  = sqrt(h P / (k Ac))   [1/m]
    base excess temp       theta0 = Tb - Tinf    [K]

Insulated tip (adiabatic tip, the primary model)::

    theta(x) = theta0 * cosh(m (L - x)) / cosh(m L)
    q_f      = sqrt(h P k Ac) * theta0 * tanh(m L)          [W]

Corrected length ``Lc = L + Ac/P`` (accounts for tip convection while
keeping the adiabatic-tip algebra; standard engineering practice)::

    q_f,c = sqrt(h P k Ac) * theta0 * tanh(m Lc)            [W]
    eta_f = tanh(m Lc) / (m Lc)        (fin efficiency, 0 < eta <= 1)
    eps_f = q_f,c / (h Ac theta0)      (fin effectiveness)

Exact convective-tip solution (used for the L -> 0 limit check, where the
fin degenerates to a flat plate convecting from its base area)::

    H = h / (m k)
    theta(x) = theta0 [cosh(m (L - x)) + H sinh(m (L - x))]
                    / [cosh(m L) + H sinh(m L)]
    q_f = sqrt(h P k Ac) theta0 [sinh(m L) + H cosh(m L)]
                    / [cosh(m L) + H sinh(m L)]              [W]

All formulas are 1-D (Bi = h t / k << 1), steady state, constant
properties, no radiation.
"""

from __future__ import annotations

import math
from typing import Union

import numpy as np

from .validation import require_finite, require_positive

__all__ = ["StraightRectangularFin"]

ArrayLike = Union[float, np.ndarray]


class StraightRectangularFin:
    """Analytic straight rectangular fin with insulated or convective tip.

    Parameters
    ----------
    length, width, thickness : float
        Fin dimensions [m]; all strictly positive.
    conductivity : float
        Fin material conductivity ``k`` [W/(m K)], strictly positive.
    h : float
        Uniform convection coefficient [W/(m2 K)], strictly positive.
    base_temperature, fluid_temperature : float
        ``T_b`` and ``T_inf`` [degC] (any finite values; ``T_b != T_inf``
        required for a meaningful fin).
    """

    def __init__(
        self,
        *,
        length: float,
        width: float,
        thickness: float,
        conductivity: float,
        h: float,
        base_temperature: float,
        fluid_temperature: float,
    ) -> None:
        self.length = require_positive("length", length)
        self.width = require_positive("width", width)
        self.thickness = require_positive("thickness", thickness)
        self.conductivity = require_positive("conductivity", conductivity)
        self.h = require_positive("h", h)
        self.base_temperature = require_finite("base_temperature", base_temperature)
        self.fluid_temperature = require_finite("fluid_temperature", fluid_temperature)

    # -- derived quantities -------------------------------------------------
    @property
    def perimeter(self) -> float:
        """Exposed perimeter ``P = 2w + 2t`` [m]."""
        return 2.0 * (self.width + self.thickness)

    @property
    def cross_section_area(self) -> float:
        """Cross-sectional (conduction) area ``Ac = w t`` [m2]."""
        return self.width * self.thickness

    @property
    def fin_parameter(self) -> float:
        """``m = sqrt(h P / (k Ac))`` [1/m]."""
        return math.sqrt(self.h * self.perimeter / (self.conductivity * self.cross_section_area))

    @property
    def base_excess_temperature(self) -> float:
        """``theta0 = T_b - T_inf`` [K]."""
        return self.base_temperature - self.fluid_temperature

    @property
    def m(self) -> float:
        """Alias of :attr:`fin_parameter` [1/m]."""
        return self.fin_parameter

    @property
    def corrected_length(self) -> float:
        """``Lc = L + Ac/P`` [m] (tip-convection corrected length)."""
        return self.length + self.cross_section_area / self.perimeter

    def _check_x(self, x: ArrayLike, *, from_tip: bool = False) -> np.ndarray:
        arr = np.asarray(x, dtype=float)
        if arr.size and (not np.all(np.isfinite(arr))):
            raise ValueError("position must be finite")
        upper = self.length if not from_tip else self.length
        if arr.size and (np.any(arr < -0.0) or np.any(arr > upper)):
            raise ValueError(f"position must lie in [0, {upper}] m")
        return arr

    # -- temperature field --------------------------------------------------
    def temperature_excess(self, x: ArrayLike) -> ArrayLike:
        """Excess temperature ``theta(x) = T(x) - T_inf`` [K], insulated tip.

        ``theta(x) = theta0 cosh(m (L - x)) / cosh(m L)`` for ``x`` in
        ``[0, L]`` measured from the base.  Scalar in, scalar out.
        """
        xx = self._check_x(x)
        m = self.m
        values = self.base_excess_temperature * np.cosh(m * (self.length - xx)) / math.cosh(m * self.length)
        return float(values) if np.ndim(x) == 0 else values

    def temperature_excess_from_tip(self, distance_from_tip: ArrayLike) -> ArrayLike:
        """Same field, parameterized by distance ``y`` from the tip [m].

        ``theta = theta0 cosh(m y) / cosh(m L)`` with ``y = L - x``; used
        for the midpoint symmetry check (two independent evaluations).
        """
        yy = self._check_x(distance_from_tip)
        m = self.m
        values = self.base_excess_temperature * np.cosh(m * yy) / math.cosh(m * self.length)
        return float(values) if np.ndim(distance_from_tip) == 0 else values

    def temperature(self, x: ArrayLike) -> ArrayLike:
        """Fin temperature ``T(x)`` [degC], insulated tip."""
        theta = self.temperature_excess(x)
        return theta + self.fluid_temperature

    def temperature_excess_convective_tip(self, x: ArrayLike) -> ArrayLike:
        """Excess temperature [K] for the exact convective-tip solution."""
        xx = self._check_x(x)
        m = self.m
        big_h = self.h / (m * self.conductivity)
        cosh_ml = math.cosh(m * self.length)
        sinh_ml = math.sinh(m * self.length)
        numerator = np.cosh(m * (self.length - xx)) + big_h * np.sinh(m * (self.length - xx))
        denominator = cosh_ml + big_h * sinh_ml
        values = self.base_excess_temperature * numerator / denominator
        return float(values) if np.ndim(x) == 0 else values

    # -- heat rates ---------------------------------------------------------
    def _m_factor(self) -> float:
        """``sqrt(h P k Ac)`` [W/K]."""
        return math.sqrt(self.h * self.perimeter * self.conductivity * self.cross_section_area)

    def heat_rate(self) -> float:
        """Fin heat rate [W], insulated tip: ``sqrt(h P k Ac) theta0 tanh(m L)``."""
        return self._m_factor() * self.base_excess_temperature * math.tanh(self.m * self.length)

    def heat_rate_corrected(self) -> float:
        """Fin heat rate [W] using the corrected length ``Lc = L + Ac/P``."""
        return self._m_factor() * self.base_excess_temperature * math.tanh(self.m * self.corrected_length)

    def heat_rate_convective_tip(self) -> float:
        """Fin heat rate [W], exact convective-tip solution."""
        m = self.m
        big_h = self.h / (m * self.conductivity)
        cosh_ml = math.cosh(m * self.length)
        sinh_ml = math.sinh(m * self.length)
        return (
            self._m_factor()
            * self.base_excess_temperature
            * (sinh_ml + big_h * cosh_ml)
            / (cosh_ml + big_h * sinh_ml)
        )

    def section_heat_rate(self, x: ArrayLike) -> ArrayLike:
        """Conductive heat rate ``q(x)`` crossing section ``x`` [W], insulated tip.

        ``q(x) = sqrt(h P k Ac) theta0 sinh(m (L - x)) / cosh(m L)``; equals
        :meth:`heat_rate` at ``x = 0`` and 0 at the insulated tip ``x = L``.
        """
        xx = self._check_x(x)
        m = self.m
        values = (
            self._m_factor()
            * self.base_excess_temperature
            * np.sinh(m * (self.length - xx))
            / math.cosh(m * self.length)
        )
        return float(values) if np.ndim(x) == 0 else values

    # -- performance metrics ------------------------------------------------
    def efficiency(self) -> float:
        """Fin efficiency ``eta_f = q_f / (h Af theta0)`` with ``Af = P Lc``.

        Reduces analytically to ``tanh(m Lc) / (m Lc)``, so ``0 < eta <= 1``.
        """
        if self.base_excess_temperature == 0.0:
            raise ValueError("efficiency is undefined when theta0 = 0")
        m_lc = self.m * self.corrected_length
        return math.tanh(m_lc) / m_lc

    def effectiveness(self) -> float:
        """Fin effectiveness ``eps_f = q_f / (h Ac theta0)`` (base-area basis)."""
        if self.base_excess_temperature == 0.0:
            raise ValueError("effectiveness is undefined when theta0 = 0")
        return self.heat_rate_corrected() / (
            self.h * self.cross_section_area * self.base_excess_temperature
        )

    # -- misc ---------------------------------------------------------------
    def __repr__(self) -> str:  # pragma: no cover - convenience only
        return (
            f"StraightRectangularFin(L={self.length} m, w={self.width} m, "
            f"t={self.thickness} m, k={self.conductivity} W/(m K), h={self.h} W/(m2 K))"
        )
