"""Composite-wall network tests: hand values, matrix vs reduction, energy balance."""

import numpy as np
import pytest

from thermal_pathways import (
    ConductionSegment,
    ContactSegment,
    ConvectionSegment,
    ParallelSegment,
    SeriesBranch,
    parallel_equivalent_resistance,
    solve_wall,
)


def recorded_wall_segments():
    """The recorded five-element wall (convection / brick / contact / insulation / convection)."""
    area = 2.0
    return [
        ConvectionSegment("indoor film", h=15.0, area=area),
        ConductionSegment("brick", thickness=0.10, conductivity=0.72, area=area),
        ContactSegment("mortar joint", resistance=0.010),
        ConductionSegment("insulation", thickness=0.05, conductivity=0.035, area=area),
        ConvectionSegment("outdoor film", h=40.0, area=area),
    ]


class TestSeriesWallHandComputed:
    def test_two_layer_exact_hand_values(self):
        """Hand-computed two-layer case (numbers stated here).

        Layer A: L = 0.02 m, k = 1.0 W/(m K), A = 1.0 m2  ->  R_A = 0.02 K/W
        Layer B: L = 0.03 m, k = 1.0 W/(m K), A = 1.0 m2  ->  R_B = 0.03 K/W
        Hot side 100 degC, cold side 20 degC:
            R_tot = 0.05 K/W,  q = 80 / 0.05 = 1600 W (exact),
            interface T = 100 - 1600 * 0.02 = 68 degC (exact).
        """
        layer_a = ConductionSegment("layer A", 0.02, 1.0, 1.0)
        layer_b = ConductionSegment("layer B", 0.03, 1.0, 1.0)
        solution = solve_wall([layer_a, layer_b], 100.0, 20.0)

        assert solution.total_resistance == pytest.approx(0.05, rel=1e-12)
        assert solution.heat_rate == pytest.approx(1600.0, rel=1e-12)
        assert solution.node_temperatures[0] == pytest.approx(100.0, abs=0.0)
        assert solution.node_temperatures[1] == pytest.approx(68.0, rel=1e-12)
        assert solution.node_temperatures[2] == pytest.approx(20.0, abs=0.0)
        for link in solution.links:
            assert link.heat_rate == pytest.approx(1600.0, rel=1e-12)
            assert link.temperature_drop == pytest.approx(
                link.resistance * 1600.0, rel=1e-12
            )


class TestMatrixAgainstDirectReduction:
    def test_node_matrix_solution_matches_voltage_divider(self):
        """Interior node temperatures from the linear solve must equal the
        direct series reduction (hand formula T_k = T_hot - q * sum(R_i, i < k))."""
        segments = recorded_wall_segments()
        solution = solve_wall(segments, 21.0, -5.0)

        resistances = [segment.resistance for segment in segments]
        q_direct = (21.0 - (-5.0)) / sum(resistances)
        assert solution.heat_rate == pytest.approx(q_direct, rel=1e-12)

        cumulative = 0.0
        expected = [21.0]
        for resistance in resistances:
            cumulative += resistance
            expected.append(21.0 - q_direct * cumulative)
        assert solution.node_temperatures == pytest.approx(
            np.array(expected), rel=1e-12
        )

    def test_single_segment_needs_no_interior_solve(self):
        layer = ConductionSegment("slab", 0.02, 1.0, 1.0)
        solution = solve_wall([layer], 100.0, 20.0)
        assert solution.heat_rate == pytest.approx(4000.0, rel=1e-12)
        assert solution.node_temperatures.tolist() == [100.0, 20.0]


class TestEnergyBalance:
    def test_wall_dissipates_exactly_what_enters(self):
        solution = solve_wall(recorded_wall_segments(), 21.0, -5.0)
        residual = solution.energy_balance_residual()
        assert residual == pytest.approx(0.0, abs=1e-10)
        assert residual < 1e-12 * solution.heat_rate
        assert solution.hot_side_supply() == pytest.approx(
            solution.cold_side_removal(), rel=1e-12
        )


class TestParallelSuperposition:
    def test_parallel_conductances_superpose_hand_case(self):
        """Hand-computed parallel case (numbers stated here).

        Path a: L = 0.05 m, k = 1.0, A = 0.1 m2 -> R_a = 0.5 K/W
        Path b: L = 0.05 m, k = 2.0, A = 0.1 m2 -> R_b = 0.25 K/W
        With 100 K across the group: q_a = 200 W, q_b = 400 W,
        total 600 W, R_eq = 1/6 K/W (exact).
        """
        path_a = SeriesBranch(
            "path a", (ConductionSegment("a", 0.05, 1.0, 0.1),)
        )
        path_b = SeriesBranch(
            "path b", (ConductionSegment("b", 0.05, 2.0, 0.1),)
        )
        group = ParallelSegment("parallel group", (path_a, path_b))
        solution = solve_wall([group], 100.0, 0.0)

        assert solution.total_resistance == pytest.approx(1.0 / 6.0, rel=1e-12)
        assert solution.heat_rate == pytest.approx(600.0, rel=1e-12)

        links = solution.heat_rates_for_group("parallel group")
        assert len(links) == 2
        by_name = {link.name: link.heat_rate for link in links}
        assert by_name["parallel group: path a"] == pytest.approx(200.0, rel=1e-12)
        assert by_name["parallel group: path b"] == pytest.approx(400.0, rel=1e-12)
        assert sum(by_name.values()) == pytest.approx(solution.heat_rate, rel=1e-12)

    def test_parallel_helper_matches_product_formula(self):
        assert parallel_equivalent_resistance([0.5, 0.25]) == pytest.approx(
            0.25 * 0.5 / (0.25 + 0.5), rel=1e-14
        )

    def test_parallel_wall_energy_balance(self):
        """Stud wall: insulation (1.78 m2) in parallel with studs (0.22 m2)."""
        insulation = SeriesBranch(
            "insulation",
            (ConductionSegment("cavity insulation", 0.09, 0.038, 1.78),),
        )
        studs = SeriesBranch(
            "timber studs", (ConductionSegment("stud", 0.09, 0.14, 0.22),)
        )
        cavity = ParallelSegment("cavity", (insulation, studs))
        solution = solve_wall([cavity], 21.0, -5.0)
        assert solution.energy_balance_residual() < 1e-12 * solution.heat_rate
        # Branch split is proportional to conductance (superposition).
        links = solution.heat_rates_for_group("cavity")
        conductances = [1.0 / link.resistance for link in links]
        expected = solution.heat_rate * np.array(conductances) / sum(conductances)
        actual = np.array([link.heat_rate for link in links])
        assert actual == pytest.approx(expected, rel=1e-12)
