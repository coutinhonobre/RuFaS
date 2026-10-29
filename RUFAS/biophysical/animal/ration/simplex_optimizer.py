"""
Ration formulation by the Simplex method (sequential linear programming).

The ration problem of ``RationOptimizer`` has linear constraints (calcium, phosphorus, dry matter intake),
ratio constraints (NDF, forage NDF and fat as percent of dry matter, linear once multiplied by the intake) and
two non-linear ones (NASEM net energy and metabolizable protein, whose coefficients depend on the whole diet).
This module linearizes the two non-linear constraints around the current diet, solves the resulting linear
program with the primal Simplex (PuLP modelling, HiGHS solver) and repeats until the diet stops moving. The
optimal tableau of the last linear program, the model (A, b, c), shadow prices and reduced costs are exported
for sensitivity analysis. The result has the same shape as the one ``scipy.optimize.minimize`` returns, so the
rest of RuFaS does not know which solver ran.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

import highspy
import numpy as np
import numpy.typing as npt
import pulp
from scipy.optimize import OptimizeResult

from RUFAS.output_manager import OutputManager
from RUFAS.units import MeasurementUnits

SIMPLEX_METHOD_NAME = "simplex"
"""Value of the ``--ration-optimizer`` command-line option that selects this solver."""

SIMPLEX_OUTPUT_DIRECTORY = Path("output") / "simplex"
"""Directory where the tableau and model files of each formulation are written."""

RATIO_CONSTRAINTS = frozenset(
    {"NDF_constraint_lower", "NDF_constraint_upper", "forage_NDF_constraint", "fat_constraint"}
)
LINEAR_CONSTRAINTS = frozenset(
    {"calcium_constraint", "phosphorus_constraint", "DMI_constraint_lower", "DMI_constraint_upper"}
)

GREATER_EQUAL, LESS_EQUAL, EQUAL = ">=", "<=", "="
PULP_SENSES = {LESS_EQUAL: pulp.LpConstraintLE, GREATER_EQUAL: pulp.LpConstraintGE, EQUAL: pulp.LpConstraintEQ}
PULP_STATUS = {
    pulp.LpStatusOptimal: "optimal",
    pulp.LpStatusInfeasible: "infeasible",
    pulp.LpStatusUnbounded: "unbounded",
    pulp.LpStatusNotSolved: "iteration_limit",
}


@dataclass
class LinearProgram:
    """``min c.x  s.t.  A x (senses) b,  lower <= x <= upper``; ``row_kinds`` says how each row was obtained."""

    c: npt.NDArray[np.float64]
    A: npt.NDArray[np.float64]
    b: npt.NDArray[np.float64]
    senses: list[str]
    row_names: list[str]
    col_names: list[str]
    lower: npt.NDArray[np.float64]
    upper: npt.NDArray[np.float64]
    row_kinds: list[str] = field(default_factory=list)


@dataclass
class SimplexResult:
    """
    Solution of a ``LinearProgram``. ``tableau`` is the optimal tableau ``B^-1 [A | S]`` (one row per basic
    variable plus the Z row, RHS in the last column) rebuilt from the basis reported by the solver; it is empty
    when no basis could be formed. Rows include the ``max_<feed>`` rows added for finite upper bounds.
    """

    status: str
    x: npt.NDArray[np.float64]
    objective: float
    tableau: npt.NDArray[np.float64]
    tableau_row_names: list[str]
    tableau_col_names: list[str]
    row_names: list[str]
    row_senses: list[str]
    row_values: npt.NDArray[np.float64]
    row_rhs: npt.NDArray[np.float64]
    shadow_prices: npt.NDArray[np.float64]
    reduced_costs: npt.NDArray[np.float64]
    basis: list[str]
    pivots: int


class LinearProgramSolver:
    """
    Solves a ``LinearProgram`` with PuLP running the primal (conventional) Simplex of HiGHS, without presolve.

    PuLP's bundled CBC is an x86_64 binary that does not run on Apple Silicon, so HiGHS (``highspy``) is used
    on every platform. ``simplex_strategy`` 4 is the primal Simplex, 1 the dual.
    """

    def __init__(self, simplex_strategy: int = 4) -> None:
        self.simplex_strategy = simplex_strategy

    def _solver(self) -> pulp.LpSolver:
        return pulp.HiGHS(msg=False, solver="simplex", simplex_strategy=self.simplex_strategy, presolve="off")

    def solve(self, lp: LinearProgram) -> SimplexResult:
        n = len(lp.c)
        c, lower = np.asarray(lp.c, dtype=float), np.asarray(lp.lower, dtype=float)
        A, b, senses, names = self._rows_with_upper_bounds(lp)
        b_shift = b - A @ lower  # x = lower + x', x' >= 0

        # PuLP sorts variables by name: zero-padded names keep the decision-vector order.
        problem = pulp.LpProblem("racao", pulp.LpMinimize)
        xs = [pulp.LpVariable(f"x{j:04d}", lowBound=0.0) for j in range(n)]
        problem += pulp.lpSum(float(c[j]) * xs[j] for j in range(n) if c[j])
        for i in range(len(b)):
            row = pulp.lpSum(float(A[i, j]) * xs[j] for j in range(n) if A[i, j])
            problem += pulp.LpConstraint(row, PULP_SENSES[senses[i]], f"r{i:04d}", float(b_shift[i]))
        with np.errstate(all="ignore"):
            status = PULP_STATUS.get(problem.solve(self._solver()), "numerical_error")
        highs = problem.solverModel
        pivots = int(highs.getInfo().simplex_iteration_count)
        if status != "optimal":
            m, nan = len(b), float("nan")
            return SimplexResult(
                status,
                lower.copy(),
                nan,
                np.zeros((0, 0)),
                [],
                list(lp.col_names) + ["RHS"],
                names,
                senses,
                np.full(m, nan),
                np.full(m, nan),
                np.full(m, nan),
                np.full(n, nan),
                [],
                pivots,
            )

        solution, basis_status = highs.getSolution(), highs.getBasis()
        x = lower + np.asarray(list(solution.col_value), dtype=float)
        shadow = np.asarray(list(solution.row_dual), dtype=float)  # dZ/db_i in each row's own sense
        reduced = np.asarray(list(solution.col_dual), dtype=float)
        basic_cols = [j for j, s in enumerate(basis_status.col_status) if s == highspy.HighsBasisStatus.kBasic]
        basic_rows = [i for i, s in enumerate(basis_status.row_status) if s == highspy.HighsBasisStatus.kBasic]
        tableau, t_rows, t_cols, basis = self._tableau(
            A, b_shift, c, senses, names, list(lp.col_names), basic_cols, basic_rows
        )
        if tableau.size:
            tableau[-1, -1] += float(c @ lower)  # the tableau is written in the shifted variables
        return SimplexResult(
            "optimal", x, float(c @ x), tableau, t_rows, t_cols, names, senses, A @ x, b, shadow, reduced, basis, pivots
        )

    @staticmethod
    def _rows_with_upper_bounds(lp: LinearProgram):
        """Appends one ``max_<feed>`` row per finite upper bound, so that bounds get their own shadow price."""
        n = len(lp.c)
        A = [np.asarray(lp.A, dtype=float).reshape(-1, n)]
        b = [np.asarray(lp.b, dtype=float)]
        senses, names = list(lp.senses), list(lp.row_names)
        for j, upper in enumerate(lp.upper):
            if np.isfinite(upper):
                A.append(np.eye(1, n, j))
                b.append(np.array([float(upper)]))
                senses.append(LESS_EQUAL)
                names.append(f"max_{lp.col_names[j]}")
        return np.vstack(A), np.concatenate(b), senses, names

    @staticmethod
    def _tableau(A, b, c, senses, row_names, col_names, basic_cols, basic_rows):
        """Optimal tableau ``B^-1 [A | S]``: one slack column per ``<=`` row, one surplus column per ``>=`` row."""
        m, n = A.shape
        extra, extra_names, extra_of_row = [], [], {}
        for i in range(m):
            if senses[i] == EQUAL:
                continue
            extra_of_row[i] = n + len(extra)
            extra.append(np.eye(1, m, i)[0] * (1.0 if senses[i] == LESS_EQUAL else -1.0))
            extra_names.append(("folga_" if senses[i] == LESS_EQUAL else "excesso_") + row_names[i])
        A_std = np.hstack([A, np.column_stack(extra)]) if extra else A
        c_std = np.concatenate([c, np.zeros(len(extra))])
        names = list(col_names) + extra_names
        basis = basic_cols + [extra_of_row[i] for i in basic_rows if i in extra_of_row]
        basis_names = [names[k] for k in basis]
        if len(basis) != m or abs(np.linalg.det(A_std[:, basis])) < 1e-12:
            return np.zeros((0, 0)), [], names + ["RHS"], basis_names
        B_inv = np.linalg.inv(A_std[:, basis])
        body, rhs = B_inv @ A_std, B_inv @ b
        z_row, z_value = c_std - c_std[basis] @ body, float(c_std[basis] @ rhs)
        tableau = np.column_stack([np.vstack([body, z_row]), np.append(rhs, z_value)])
        return tableau, basis_names + ["Z"], names + ["RHS"], basis_names


@dataclass
class SimplexRationReport:
    """What one formulation produced, for export: last LP and its solution, constraint values, iterations."""

    lp: LinearProgram | None
    result: SimplexResult | None
    constraint_values: dict[str, float]
    slp_iterations: int
    status: str
    message: str


class SimplexRationSolver:
    """
    Sequential linear programming for the ration problem.

    Each iteration linearizes the constraints around the current diet ``x_k`` (exactly for linear and ratio
    constraints, by the tangent plane for the others), solves the linear program inside a trust region around
    ``x_k`` and accepts the new diet when the merit (cost plus a penalty on constraint violation) does not get
    worse; otherwise the trust region shrinks. It stops when the diet stops moving. A final pass re-solves the
    program at the final diet with the real feed bounds, so the exported tableau belongs to that diet, and
    tightens any constraint still violated by rounding or curvature.
    """

    def __init__(self, max_iterations: int = 60, step_tolerance: float = 1e-4, penalty_weight: float = 1_000.0) -> None:
        self.max_iterations = max_iterations
        self.step_tolerance = step_tolerance
        self.penalty_weight = penalty_weight
        self.simplex = LinearProgramSolver()

    def solve(
        self,
        objective_coefficients: Sequence[float],
        constraints: Sequence[dict[str, Any]],
        ration_config: Any,
        bounds: Sequence[tuple[float, float]],
        initial_x: Sequence[float],
        variable_names: Sequence[str],
    ) -> tuple[OptimizeResult, SimplexRationReport]:
        """Returns a scipy-like result (``x``, ``success``, ``fun``, ``nit``, ``message``) and the export report."""
        c = np.asarray(objective_coefficients, dtype=float)
        lower = np.array([max(0.0, float(lo)) for lo, _ in bounds])
        upper = np.array([float(hi) if hi is not None else np.inf for _, hi in bounds])
        names = [str(v) for v in variable_names]
        functions = [(k["fun"].__name__, k["fun"]) for k in constraints]
        scales = self._scales(functions, ration_config, len(c))

        def values_at(x: npt.NDArray[np.float64]) -> dict[str, float]:
            return {name: float(fun(np.asarray(x, dtype=float), ration_config)) for name, fun in functions}

        def merit(x: npt.NDArray[np.float64]) -> float:
            violation = max([0.0] + [-v / scales[k] for k, v in values_at(x).items()])
            return float(c @ x) + self.penalty_weight * violation

        x = np.clip(np.asarray(initial_x, dtype=float), lower, np.minimum(upper, 1e12))
        intake = max(1.0, float(x.sum()))
        tol, radius, min_radius = self.step_tolerance * intake, intake, 0.1 * self.step_tolerance * intake
        margins = np.zeros(len(functions))
        last: tuple[LinearProgram, SimplexResult] | None = None
        merit_x, failure, converged, iterations = merit(x), None, False, 0

        for iterations in range(1, self.max_iterations + 1):
            lp = self._build_lp(c, functions, ration_config, x, lower, upper, radius, margins, names)
            result = self.simplex.solve(lp)
            if result.status != "optimal":
                if radius < 4.0 * intake:  # the trust region may exclude the feasible set: widen it
                    radius = min(2.0 * radius, 4.0 * intake)
                    continue
                last, failure = (lp, result), result.status
                break
            step, merit_new = float(np.abs(result.x - x).max()), merit(result.x)
            if merit_new > merit_x + 1e-9 * max(1.0, abs(merit_x)):  # linearization over-promised: reject
                if radius <= min_radius:
                    converged = True
                    break
                radius = max(radius / 2.0, min_radius)
                continue
            last, stalled = (lp, result), merit_x - merit_new <= 1e-8 * max(1.0, abs(merit_x))
            x, merit_x = result.x, merit_new
            at_border = step >= 0.9 * radius
            if (step <= tol and not at_border) or (stalled and step <= 10.0 * tol):
                converged = True
                break
            if at_border:
                radius = min(2.0 * radius, 4.0 * intake)

        if failure is None:
            x, last = self._finalize(
                c, functions, ration_config, x, lower, upper, margins, names, values_at, merit, last
            )
        values = values_at(x)
        success = failure is None and all(v >= -1e-9 * scales[k] for k, v in values.items())
        if success:
            status = "optimal"
            message = "Sequential Simplex converged." if converged else "Iteration limit reached; ration is feasible."
        else:
            status = failure or "infeasible"
            message = f"Sequential Simplex failed: {status}."
        out = OptimizeResult(
            x=x, success=success, fun=float(c @ x), nit=iterations, status=int(not success), message=message
        )
        report = SimplexRationReport(
            last[0] if last else None, last[1] if last else None, values, iterations, status, message
        )
        return out, report

    @staticmethod
    def _scales(functions: Sequence[tuple[str, Callable[..., float]]], ration_config: Any, n: int) -> dict[str, float]:
        """Magnitude of each constraint (its value at the empty ration), to judge violations relatively."""
        scales = {}
        with np.errstate(all="ignore"):
            for name, fun in functions:
                try:
                    value = float(fun(np.zeros(n), ration_config))
                except (ZeroDivisionError, ValueError, FloatingPointError):
                    value = float("nan")
                scales[name] = max(1.0, abs(value)) if np.isfinite(value) else 1.0
        return scales

    def _build_lp(self, c, functions, ration_config, x_k, lower, upper, radius, margins, names) -> LinearProgram:
        """Linearizes every constraint around ``x_k`` (all rows become ``a.x >= b``) inside the trust region."""
        rows, rhs, kinds = [], [], []
        for (name, fun), margin in zip(functions, margins):
            a, b, kind = self._linearize(name, fun, ration_config, x_k)
            rows.append(a)
            rhs.append(b + margin)
            kinds.append(kind)
        lp_lower = np.maximum(lower, x_k - radius)
        lp_upper = np.maximum(np.minimum(upper, x_k + radius), lp_lower)
        row_names = [name for name, _ in functions]
        return LinearProgram(
            c,
            np.vstack(rows),
            np.asarray(rhs),
            [GREATER_EQUAL] * len(rows),
            row_names,
            names,
            lp_lower,
            lp_upper,
            kinds,
        )

    @staticmethod
    def _linearize(name: str, fun: Callable[..., float], ration_config: Any, x_k: npt.NDArray[np.float64]):
        """Returns ``(a, b, kind)`` with ``fun(x) >= 0`` represented by ``a.x >= b``."""
        n = len(x_k)

        def g(z: npt.NDArray[np.float64]) -> float:
            return float(fun(np.asarray(z, dtype=float), ration_config))

        unit = np.eye(n)
        if name in RATIO_CONSTRAINTS and x_k.sum() > 0:  # (sum a_j x_j)/(sum x_j) - k >= 0  <=>  sum (a_j - k) x_j >= 0
            a = np.array([g(unit[j]) for j in range(n)])
            if np.isclose(g(2.0 * x_k), g(x_k), rtol=1e-7, atol=1e-9):
                return a, 0.0, "razao"
        if name in LINEAR_CONSTRAINTS:  # affine: recover the coefficients exactly
            g0 = g(np.zeros(n))
            a = np.array([g(unit[j]) - g0 for j in range(n)])
            if np.isclose(g(x_k), a @ x_k + g0, rtol=1e-7, atol=1e-9):
                return a, -g0, "linear"
        gk, grad = g(x_k), np.zeros(n)  # tangent plane by finite differences
        for j in range(n):
            h = 1e-5 * max(1.0, abs(float(x_k[j])))
            fwd, bwd = x_k.copy(), x_k.copy()
            fwd[j] += h
            bwd[j] -= h
            grad[j] = (g(fwd) - g(bwd)) / (2.0 * h) if bwd[j] >= 0.0 else (g(fwd) - gk) / h
        return grad, float(grad @ x_k - gk), "linearizada"

    def _finalize(self, c, functions, ration_config, x, lower, upper, margins, names, values_at, merit, last):
        """Re-solves at the final diet with the real bounds (for the exported tableau), then repairs violations."""

        def worst(z: npt.NDArray[np.float64]) -> float:
            return max([0.0] + [-v for v in values_at(z).values()])

        lp = self._build_lp(c, functions, ration_config, x, lower, upper, np.inf, margins, names)
        result = self.simplex.solve(lp)
        if result.status == "optimal":
            if merit(result.x) <= merit(x) + 1e-9 * max(1.0, abs(merit(x))) and worst(result.x) <= worst(x):
                x = result.x
            last = (lp, result)
        for _ in range(8):  # tighten violated rows and re-solve near x until every original constraint holds
            values = values_at(x)
            violated = [i for i, (name, _) in enumerate(functions) if values[name] < 0.0]
            if not violated:
                break
            for i in violated:
                margins[i] = max(2.0 * margins[i], 2.0 * -values[functions[i][0]] + 1e-7)
            small = 1e-3 * max(1.0, float(x.sum()))
            repair = self.simplex.solve(
                self._build_lp(c, functions, ration_config, x, lower, upper, small, margins, names)
            )
            if repair.status != "optimal":
                break
            if worst(repair.x) < worst(x):
                x = repair.x
        return x, last


class SimplexReportExporter:
    """Writes the tableau and the model of a formulation to CSV files and the summary to the ``OutputManager``."""

    _attempts: dict[tuple[str, int], int] = {}

    @classmethod
    def export(cls, pen_label: str, report: SimplexRationReport, feed_ids: Sequence[Any]) -> Path | None:
        om = OutputManager()
        day = int(om.time.simulation_day) if om.time is not None else 0
        cls._attempts[(pen_label, day)] = attempt = cls._attempts.get((pen_label, day), 0) + 1
        info = {
            "class": cls.__name__,
            "function": cls.export.__name__,
            "simulation_day": day,
            "is_daily_variable": True,
        }
        result = report.result
        summary = {
            "status": report.status,
            "slp_iterations": report.slp_iterations,
            "attempt": attempt,
            "objective": result.objective if result else float("nan"),
            "pivots": result.pivots if result else 0,
        }
        units = {k: (MeasurementUnits.DOLLARS if k == "objective" else MeasurementUnits.UNITLESS) for k in summary}
        om.add_variable(f"simplex_for_{pen_label}", summary, {**info, "units": units})
        if report.lp is None or result is None or result.status != "optimal":
            return None
        shadow = {name: float(v) for name, v in zip(result.row_names, result.shadow_prices)}
        reduced = {str(f): float(v) for f, v in zip(feed_ids, result.reduced_costs)}
        for name, values in (
            (f"simplex_shadow_price_for_{pen_label}", shadow),
            (f"simplex_reduced_cost_for_{pen_label}", reduced),
        ):
            om.add_variable(name, values, {**info, "units": {k: MeasurementUnits.DOLLARS for k in values}})

        prefix = getattr(om, "_OutputManager__metadata_prefix", "") or "simulation"  # private in OutputManager
        directory = SIMPLEX_OUTPUT_DIRECTORY / prefix
        directory.mkdir(parents=True, exist_ok=True)
        base = f"{pen_label}_dia{day:05d}_tent{attempt}"
        cls._write_tableau(directory / f"tableau_{base}.csv", result)
        cls._write_model(directory / f"modelo_{base}.csv", report, feed_ids)
        return directory / f"tableau_{base}.csv"

    @staticmethod
    def _write_tableau(path: Path, result: SimplexResult) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["base"] + result.tableau_col_names)
            for name, row in zip(result.tableau_row_names, result.tableau):
                writer.writerow([name] + [f"{v:.10g}" for v in row])

    @staticmethod
    def _write_model(path: Path, report: SimplexRationReport, feed_ids: Sequence[Any]) -> None:
        lp, result = report.lp, report.result
        assert lp is not None and result is not None
        feeds = [str(f) for f in feed_ids]

        def fmt(values: Sequence[float]) -> list[str]:
            return [f"{v:.10g}" for v in values]

        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["restricao", "tipo"] + feeds + ["sinal", "b", "valor", "folga", "preco_sombra"])
            for i, name in enumerate(result.row_names):
                if i < lp.A.shape[0]:
                    coefficients, kind = lp.A[i], lp.row_kinds[i]
                else:  # upper-bound row added by the solver
                    bounded = name.removeprefix("max_")
                    coefficients, kind = np.array([1.0 if f == bounded else 0.0 for f in feeds]), "limite_superior"
                slack = result.row_values[i] - result.row_rhs[i]
                tail = fmt([result.row_rhs[i], result.row_values[i], slack, result.shadow_prices[i]])
                writer.writerow([name, kind] + fmt(coefficients) + [result.row_senses[i]] + tail)
            writer.writerow([])
            writer.writerow(["variavel", ""] + feeds)
            rows = (
                ("custo_c", lp.c),
                ("x_otimo", result.x),
                ("custo_reduzido", result.reduced_costs),
                ("limite_inferior", lp.lower),
                ("limite_superior", lp.upper),
            )
            for label, values in rows:
                writer.writerow([label, ""] + fmt(values))
            writer.writerow([])
            writer.writerow(["resumo", "valor"])
            writer.writerow(["Z", f"{result.objective:.10g}"])
            writer.writerow(["status", report.status])
            writer.writerow(["iteracoes_linearizacao", report.slp_iterations])
            writer.writerow(["pivos_simplex", result.pivots])
            writer.writerow(["base_otima", " ".join(result.basis)])
            writer.writerow([])
            writer.writerow(["restricao_original", "valor_em_x_otimo"])
            for name, value in report.constraint_values.items():
                writer.writerow([name, f"{value:.10g}"])
