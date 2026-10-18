"""Input validation helpers shared across :mod:`thermal_pathways`.

Every physical quantity that must be strictly positive (dimensions,
conductivities, heat-transfer coefficients, areas, resistances) is checked
with :func:`require_positive`.  Temperatures may be negative (Celsius), so
they are only checked to be finite real numbers via :func:`require_finite`.
NaN and infinity are rejected everywhere.
"""

from __future__ import annotations

import math
from typing import Any

__all__ = ["as_float", "require_positive", "require_finite"]


def as_float(name: str, value: Any) -> float:
    """Return ``value`` as a float after checking it is a real, finite number.

    Raises
    ------
    TypeError
        If ``value`` is not a real number (``bool`` is rejected explicitly).
    ValueError
        If ``value`` is NaN or infinite.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real number, got {type(value).__name__}")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return value


def require_positive(name: str, value: Any) -> float:
    """Validate that ``value`` is a finite, strictly positive float.

    Used for lengths, widths, thicknesses, areas, conductivities ``k``,
    heat-transfer coefficients ``h`` and resistances ``R``.
    """
    value = as_float(name, value)
    if value <= 0.0:
        raise ValueError(f"{name} must be strictly positive, got {value}")
    return value


def require_finite(name: str, value: Any) -> float:
    """Validate that ``value`` is a finite float (sign is allowed)."""
    return as_float(name, value)
