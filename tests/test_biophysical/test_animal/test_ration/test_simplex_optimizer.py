from typing import Any

import numpy as np
import pytest
from scipy.optimize import linprog

from RUFAS.biophysical.animal.ration.simplex_optimizer import (
    GREATER_EQUAL,
    LESS_EQUAL,
    LinearProgram,
    SimplexRationSolver,
    LinearProgramSolver,
)


def _lp(
    c: list[float], A: list[list[float]], b: list[float], senses: list[str], lower=None, upper=None
) -> LinearProgram:
    n = len(c)
    return LinearProgram(
        c=np.array(c, dtype=float),
        A=np.array(A, dtype=float),
        b=np.array(b, dtype=float),
        senses=senses,
        row_names=[f"r{i}" for i in range(len(b))],
        col_names=[f"x{j}" for j in range(n)],
        lower=np.zeros(n) if lower is None else np.array(lower, dtype=float),
        upper=np.full(n, np.inf) if upper is None else np.array(upper, dtype=float),
    )


def test_linear_program_solver_textbook_problem() -> None:
    # min 2x + 3y  s.t.  x + y >= 4,  x <= 3,  y >= 0.5  ->  x = 3, y = 1, Z = 9
    lp = _lp([2, 3], [[1, 1], [1, 0], [0, 1]], [4, 3, 0.5], [GREATER_EQUAL, LESS_EQUAL, GREATER_EQUAL])
    result = LinearProgramSolver().solve(lp)
    assert result.status == "optimal"
    assert np.allclose(result.x, [3, 1])
    assert np.isclose(result.objective, 9)
    # Rebuilt optimal tableau: one row per basic variable plus the Z row, RHS as last column.
    assert result.tableau.shape == (4, len(result.tableau_col_names))
    assert result.tableau_row_names[-1] == "Z"
    assert np.isclose(result.tableau[-1, -1], 9)
    assert set(result.basis) == {"x0", "x1", "excesso_r2"}
    # Reduced costs are non-negative at a minimum and zero for basic variables.
    assert np.all(result.reduced_costs >= -1e-9)
    assert np.allclose(result.reduced_costs, [0, 0])
    # One more unit of demand is supplied by y (x is capped), so the shadow price is 3 $/unit;
    # relaxing the cap on x saves 1 $/unit.
    assert np.isclose(result.shadow_prices[0], 3)
    assert np.isclose(result.shadow_prices[1], -1)


def test_linear_program_solver_detects_infeasibility() -> None:
    lp = _lp([1, 1], [[1, 1], [1, 1]], [5, 3], [GREATER_EQUAL, LESS_EQUAL])
    result = LinearProgramSolver().solve(lp)
    assert result.status == "infeasible"


def test_linear_program_solver_detects_unboundedness() -> None:
    lp = _lp([-1, 0], [[0, 1]], [1], [LESS_EQUAL])
    result = LinearProgramSolver().solve(lp)
    assert result.status == "unbounded"


def test_linear_program_solver_handles_lower_and_upper_bounds() -> None:
    # min x + y  s.t.  x + y >= 2,  0.5 <= x <= 0.75,  1 <= y <= 5  -> x = 0.75 or any split; Z = 2
    lp = _lp([1, 1], [[1, 1]], [2], [GREATER_EQUAL], lower=[0.5, 1.0], upper=[0.75, 5.0])
    result = LinearProgramSolver().solve(lp)
    assert result.status == "optimal"
    assert np.isclose(result.objective, 2)
    assert 0.5 - 1e-9 <= result.x[0] <= 0.75 + 1e-9
    assert result.x[1] >= 1 - 1e-9
    assert "max_x0" in result.row_names and "max_x1" in result.row_names


def test_linear_program_solver_matches_scipy_on_campo_dourado_instance() -> None:
    # Instance of MODELAGEM_SIMPLEX_CAMPO_DOURADO.md, section 5 (6 feeds, 11 constraints).
    c = [0.16, 0.18, 0.19, 0.41, 0.36, 0.80]
    A = [
        [1.75, 1.54, 2.06, 2.36, 1.87, 0],
        [77.8, 92.1, 96.8, 236.5, 102.6, 4.7],
        [77.8, 92.1, 96.8, 236.5, 102.6, 4.7],
        [0.71, 3.76, 0.21, 2.40, 1.04, 146.6],
        [1.48, 2.25, 2.16, 5.17, 4.35, 18.56],
        [15.93, 18.20, -15.24, -13.94, 25.58, -24.86],
        [4.07, 1.80, 35.24, 33.94, -5.58, 44.86],
        [25.93, 28.20, -15, -15, -15, -15],
        [4.08, 4.13, 3.15, 5.19, -11.62, 7.00],
        [1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1],
    ]
    b = [45.5, 2575, 3863, 65.2, 61.7, 0, 0, 0, 0, 17.0, 31.7]
    senses = [GREATER_EQUAL] * 2 + [LESS_EQUAL] + [GREATER_EQUAL] * 7 + [LESS_EQUAL]
    lp = _lp(c, A, b, senses)
    result = LinearProgramSolver().solve(lp)
    assert result.status == "optimal"

    A_ub, b_ub = [], []
    for row, rhs, sense in zip(A, b, senses):
        if sense == LESS_EQUAL:
            A_ub.append(row)
            b_ub.append(rhs)
        else:
            A_ub.append([-v for v in row])
            b_ub.append(-rhs)
    reference = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(0, None)] * 6, method="highs")
    assert reference.success
    assert np.isclose(result.objective, reference.fun, atol=1e-6)
    values = np.array(A) @ result.x
    for value, rhs, sense in zip(values, b, senses):
        if sense == LESS_EQUAL:
            assert value <= rhs + 1e-7
        else:
            assert value >= rhs - 1e-7
    # Duals from the tableau agree with scipy's (sign convention: dZ/db).
    scipy_duals = -reference.ineqlin.marginals
    scipy_duals = [d if s == LESS_EQUAL else -d for d, s in zip(scipy_duals, senses)]
    assert np.allclose(result.shadow_prices[: len(b)], [-d for d in scipy_duals], atol=1e-6) or np.allclose(
        result.shadow_prices[: len(b)], scipy_duals, atol=1e-6
    )


class _Config:
    """Stand-in for RationConfig: only the attributes the fake constraints read."""

    price_list = [1.0, 2.0, 3.0]
    requirement = 10.0


def calcium_constraint(x: np.ndarray, config: Any) -> float:
    return float(np.array([1.0, 2.0, 0.5]) @ x - config.requirement)


def NDF_constraint_lower(x: np.ndarray, config: Any) -> float:
    total = float(np.sum(x))
    return float(np.array([40.0, 10.0, 60.0]) @ x / total - 25.0) if total else -1.0


def DMI_constraint_upper(x: np.ndarray, config: Any) -> float:
    return float(12.0 - np.sum(x))


def energy_constraint(x: np.ndarray, config: Any) -> float:
    # Non-linear: energy density decreases with total intake, like the NASEM discount.
    total = float(np.sum(x))
    density = np.array([2.0, 1.5, 2.5]) * (1.0 - 0.02 * total)
    return float(density @ x - 12.0)


def test_sequential_simplex_solves_mixed_constraint_problem() -> None:
    constraints = [
        {"type": "ineq", "fun": calcium_constraint},
        {"type": "ineq", "fun": NDF_constraint_lower},
        {"type": "ineq", "fun": DMI_constraint_upper},
        {"type": "ineq", "fun": energy_constraint},
    ]
    config = _Config()
    result, report = SimplexRationSolver().solve(
        objective_coefficients=config.price_list,
        constraints=constraints,
        ration_config=config,
        bounds=[(0.0, 100.0)] * 3,
        initial_x=[3.0, 3.0, 3.0],
        variable_names=["a", "b", "c"],
    )
    assert result.success, result.message
    for constraint in constraints:
        assert constraint["fun"](result.x, config) >= 0.0
    assert report.lp is not None and report.result is not None
    assert report.lp.row_kinds == ["linear", "razao", "linear", "linearizada"]
    assert report.result.status == "optimal"
    # The returned cost is the cheapest found; it must not exceed the cost of the starting point.
    assert result.fun <= float(np.array(config.price_list) @ np.array([3.0, 3.0, 3.0])) + 1e-9


def test_sequential_simplex_reports_infeasible_problem() -> None:
    constraints = [
        {"type": "ineq", "fun": calcium_constraint},
        {"type": "ineq", "fun": lambda x, config: float(1.0 - np.sum(x))},  # total <= 1 cannot reach calcium 10
    ]
    config = _Config()
    result, report = SimplexRationSolver().solve(
        objective_coefficients=config.price_list,
        constraints=constraints,
        ration_config=config,
        bounds=[(0.0, 100.0)] * 3,
        initial_x=[1.0, 1.0, 1.0],
        variable_names=["a", "b", "c"],
    )
    assert not result.success
    assert report.status in ("infeasible", "iteration_limit")
    assert result.x.shape == (3,)
