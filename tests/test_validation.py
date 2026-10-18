"""Input-validation tests (positive/finite checks, NaN/inf rejection)."""

import math

import pytest

from thermal_pathways import (
    ConductionSegment,
    ContactSegment,
    ConvectionSegment,
    ParallelSegment,
    SeriesBranch,
    StraightRectangularFin,
    require_finite,
    require_positive,
    solve_wall,
)


class TestScalarValidation:
    def test_positive_accepts_valid(self):
        assert require_positive("x", 2) == 2.0
        assert require_positive("x", 1e-9) == 1e-9

    @pytest.mark.parametrize("bad", [0.0, -1.5, -1e-9, math.inf, -math.inf, math.nan])
    def test_positive_rejects(self, bad):
        with pytest.raises(ValueError):
            require_positive("x", bad)

    def test_positive_rejects_non_numeric(self):
        with pytest.raises(TypeError):
            require_positive("x", "3")
        with pytest.raises(TypeError):
            require_positive("x", True)

    def test_finite_allows_negative_temperatures(self):
        assert require_finite("T", -40.0) == -40.0

    @pytest.mark.parametrize("bad", [math.nan, math.inf])
    def test_finite_rejects_non_finite(self, bad):
        with pytest.raises(ValueError):
            require_finite("T", bad)


class TestSegmentValidation:
    def test_conduction_rejects_non_positive(self):
        with pytest.raises(ValueError):
            ConductionSegment("bad", thickness=0.0, conductivity=1.0, area=1.0)
        with pytest.raises(ValueError):
            ConductionSegment("bad", thickness=0.1, conductivity=-1.0, area=1.0)
        with pytest.raises(ValueError):
            ConductionSegment("bad", thickness=0.1, conductivity=1.0, area=math.nan)

    def test_convection_rejects_non_positive_h(self):
        with pytest.raises(ValueError):
            ConvectionSegment("bad", h=0.0, area=1.0)

    def test_contact_rejects_non_positive_resistance(self):
        with pytest.raises(ValueError):
            ContactSegment("bad", resistance=0.0)

    def test_parallel_needs_two_branches(self):
        branch = SeriesBranch("only", (ConductionSegment("a", 0.02, 1.0, 1.0),))
        with pytest.raises(ValueError):
            ParallelSegment("solo", (branch,))

    def test_series_branch_rejects_empty_and_nested_parallel(self):
        with pytest.raises(ValueError):
            SeriesBranch("empty", ())
        nested = ParallelSegment(
            "nested",
            (
                SeriesBranch("a", (ConductionSegment("a", 0.02, 1.0, 1.0),)),
                SeriesBranch("b", (ConductionSegment("b", 0.02, 1.0, 1.0),)),
            ),
        )
        with pytest.raises(ValueError):
            SeriesBranch("outer", (nested,))


class TestFinValidation:
    def test_fin_rejects_non_positive_dimensions(self):
        kwargs = dict(
            length=0.1,
            width=0.05,
            thickness=0.003,
            conductivity=180.0,
            h=60.0,
            base_temperature=95.0,
            fluid_temperature=25.0,
        )
        for key in ("length", "width", "thickness", "conductivity", "h"):
            with pytest.raises(ValueError):
                StraightRectangularFin(**{**kwargs, key: 0.0})
        with pytest.raises(ValueError):
            StraightRectangularFin(**{**kwargs, "conductivity": math.inf})
        with pytest.raises(ValueError):
            StraightRectangularFin(**{**kwargs, "base_temperature": math.nan})

    def test_fin_rejects_position_outside_domain(self):
        fin = StraightRectangularFin(
            length=0.1,
            width=0.05,
            thickness=0.003,
            conductivity=180.0,
            h=60.0,
            base_temperature=95.0,
            fluid_temperature=25.0,
        )
        with pytest.raises(ValueError):
            fin.temperature_excess(0.11)
        with pytest.raises(ValueError):
            fin.section_heat_rate(-1e-6)
        with pytest.raises(ValueError):
            fin.section_heat_rate(math.nan)


class TestSolveWallValidation:
    def test_solve_wall_requires_segments(self):
        with pytest.raises(ValueError):
            solve_wall([], 100.0, 20.0)

    def test_solve_wall_rejects_bad_temperature_order(self):
        layer = ConductionSegment("slab", 0.02, 1.0, 1.0)
        with pytest.raises(ValueError):
            solve_wall([layer], 20.0, 100.0)
        with pytest.raises(ValueError):
            solve_wall([layer], math.nan, 20.0)
