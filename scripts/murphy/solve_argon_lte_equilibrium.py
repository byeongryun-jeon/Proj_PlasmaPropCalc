#!/usr/bin/env python3
"""Solve LTE equilibrium number densities for Ar plasma species.

Unknowns:
    n_Ar, n_Ar+, n_Ar2+, n_Ar3+, n_Ar4+, n_e

Equations:
    - 4 Saha-equilibrium equations (mu_z = mu_{z+1} + mu_e)
    - Quasi-neutrality
    - Equation of state with Debye-Huckel pressure correction
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Sequence

try:
    from scipy.optimize import root as scipy_root  # type: ignore
except ModuleNotFoundError:
    scipy_root = None


# --- Physical constants (SI) ---
K_B = 1.380649e-23  # J/K
H_PLANCK = 6.62607015e-34  # J s
EPS0 = 8.8541878128e-12  # F/m
E_CHARGE = 1.602176634e-19  # C
EV_TO_J = E_CHARGE
ATM_TO_PA = 101325.0

# --- Species definitions ---
SPECIES = ["Ar", "Ar+", "Ar2+", "Ar3+", "Ar4+", "e-"]
CHARGES = [0, 1, 2, 3, 4, -1]
ABS_CHARGES = [0, 1, 2, 3, 4, 1]

MASS_AR = 6.6335209e-26  # kg
MASS_E = 9.1093837e-31  # kg

MIN_LN_N = -5000.0
MAX_LN_N = 700.0
MIN_POSITIVE_N = 1e-300


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Solve Ar LTE equilibrium composition over T-P grid.")
    parser.add_argument(
        "--partition-json",
        type=Path,
        default=root / "data" / "processed" / "partition_functions" / "partition_functions_Ar_system.json",
        help="Path to partition-function JSON generated in Phase 1.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "data" / "processed" / "equilibrium",
        help="Directory for equilibrium solver outputs.",
    )
    parser.add_argument("--t-min", type=float, default=300.0)
    parser.add_argument("--t-max", type=float, default=30000.0)
    parser.add_argument("--t-step", type=float, default=200.0)
    parser.add_argument(
        "--pressures-atm",
        type=float,
        nargs="+",
        default=[0.1, 1.0, 4.0],
        help="Pressure cases in atm.",
    )
    parser.add_argument("--tol", type=float, default=1e-8, help="Residual tolerance for solver.")
    parser.add_argument("--max-iter", type=int, default=80, help="Max Newton iterations.")
    parser.add_argument(
        "--exclude-t-max",
        action="store_true",
        help="Do not force include t-max when off-grid.",
    )
    return parser.parse_args()


def build_temperature_grid(t_min: float, t_max: float, t_step: float, include_t_max: bool) -> list[float]:
    if t_step <= 0:
        raise ValueError("t-step must be > 0")
    if t_min <= 0:
        raise ValueError("t-min must be > 0")
    if t_max < t_min:
        raise ValueError("t-max must be >= t-min")
    out: list[float] = []
    t = t_min
    while t <= t_max + 1e-12:
        out.append(float(t))
        t += t_step
    if include_t_max and (not out or abs(out[-1] - t_max) > 1e-9):
        out.append(float(t_max))
    return out


def interp_linear(x_grid: Sequence[float], y_grid: Sequence[float], x: float) -> float:
    if x <= x_grid[0]:
        return float(y_grid[0])
    if x >= x_grid[-1]:
        return float(y_grid[-1])
    lo, hi = 0, len(x_grid) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x_grid[mid] <= x:
            lo = mid
        else:
            hi = mid
    x0, x1 = x_grid[lo], x_grid[hi]
    y0, y1 = y_grid[lo], y_grid[hi]
    w = (x - x0) / (x1 - x0)
    return float(y0 + w * (y1 - y0))


def thermal_de_broglie_wavelength(mass_kg: float, temperature_k: float) -> float:
    return H_PLANCK / math.sqrt(2.0 * math.pi * mass_kg * K_B * temperature_k)


def debye_huckel_pressure_correction(n: Sequence[float], temperature_k: float) -> float:
    # n order: Ar, Ar+, Ar2+, Ar3+, Ar4+, e-
    sigma = 0.0
    for ni, zi_abs in zip(n, ABS_CHARGES):
        sigma += (zi_abs * zi_abs) * ni
    if sigma <= 0.0:
        return 0.0
    lambda_d = math.sqrt(EPS0 * K_B * temperature_k / (E_CHARGE * E_CHARGE * sigma))
    if lambda_d <= 0.0:
        return 0.0
    kappa = 1.0 / lambda_d
    return K_B * temperature_k * (kappa**3) / (24.0 * math.pi)


def l2_norm(vec: Sequence[float]) -> float:
    return math.sqrt(sum(v * v for v in vec))


def exp_clamped(x: float) -> float:
    # For very negative log-density, return 0 in number-density space while
    # keeping the log-state itself finite for Saha relations.
    if x < -745.0:
        return 0.0
    if x > MAX_LN_N:
        return math.exp(MAX_LN_N)
    return math.exp(x)


def ln_clamped_from_n(n: float) -> float:
    return math.log(max(n, MIN_POSITIVE_N))


def lin_solve_gauss(
    a: list[list[float]], b: list[float], pivot_tol: float = 1e-40
) -> list[float]:
    n = len(a)
    m = [row[:] + [rhs] for row, rhs in zip(a, b)]

    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(m[i][k]))
        if abs(m[pivot][k]) < pivot_tol:
            raise RuntimeError("Singular Jacobian in Newton solve.")
        if pivot != k:
            m[k], m[pivot] = m[pivot], m[k]

        piv = m[k][k]
        for j in range(k, n + 1):
            m[k][j] /= piv
        for i in range(k + 1, n):
            fac = m[i][k]
            if fac == 0.0:
                continue
            for j in range(k, n + 1):
                m[i][j] -= fac * m[k][j]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = m[i][n]
        for j in range(i + 1, n):
            s -= m[i][j] * x[j]
        x[i] = s
    return x


def solve_linear_regularized(a: list[list[float]], b: list[float]) -> list[float]:
    reg_levels = [0.0, 1e-14, 1e-12, 1e-10, 1e-8]
    for reg in reg_levels:
        trial = [row[:] for row in a]
        if reg > 0.0:
            for i in range(len(trial)):
                trial[i][i] += reg
        try:
            return lin_solve_gauss(trial, b, pivot_tol=1e-40)
        except RuntimeError:
            continue
    raise RuntimeError("Singular Jacobian in Newton solve (even after regularization).")


class ArgonLTEEquilibriumSolver:
    def __init__(self, partition_json_path: Path) -> None:
        payload = json.loads(partition_json_path.read_text(encoding="utf-8"))
        self.temperature_grid = [float(t) for t in payload["temperature_K"]]
        q_payload = payload["Q"]
        self.q_grid = {sp: [float(v) for v in q_payload[sp]] for sp in SPECIES}

        eion = payload.get("cumulative_ionization_eV_by_species")
        if not eion:
            raise ValueError(
                "Partition JSON does not contain cumulative ionization energies. "
                "Rebuild partition data with updated parser."
            )
        self.eion_j = {sp: float(eion[sp]) * EV_TO_J for sp in SPECIES}

        self.masses = {
            "Ar": MASS_AR,
            "Ar+": MASS_AR - 1.0 * MASS_E,
            "Ar2+": MASS_AR - 2.0 * MASS_E,
            "Ar3+": MASS_AR - 3.0 * MASS_E,
            "Ar4+": MASS_AR - 4.0 * MASS_E,
            "e-": MASS_E,
        }

    def q_int(self, species: str, temperature_k: float) -> float:
        return interp_linear(self.temperature_grid, self.q_grid[species], temperature_k)

    def saha_log_constants(self, temperature_k: float) -> list[float]:
        """Return log(K_z) for reactions Ar^(z+) <-> Ar^((z+1)+) + e."""
        kt = K_B * temperature_k
        lam_e = thermal_de_broglie_wavelength(self.masses["e-"], temperature_k)
        q_e = max(self.q_int("e-", temperature_k), MIN_POSITIVE_N)

        out: list[float] = []
        for z in range(4):
            sp_z = SPECIES[z]
            sp_zp1 = SPECIES[z + 1]
            lam_z = thermal_de_broglie_wavelength(self.masses[sp_z], temperature_k)
            lam_zp1 = thermal_de_broglie_wavelength(self.masses[sp_zp1], temperature_k)
            q_z = max(self.q_int(sp_z, temperature_k), MIN_POSITIVE_N)
            q_zp1 = max(self.q_int(sp_zp1, temperature_k), MIN_POSITIVE_N)
            delta_e = self.eion_j[sp_zp1] + self.eion_j["e-"] - self.eion_j[sp_z]
            log_k = (
                math.log(q_zp1)
                + math.log(q_e)
                - math.log(q_z)
                + 3.0 * (math.log(lam_z) - math.log(lam_zp1) - math.log(lam_e))
                - delta_e / kt
            )
            out.append(log_k)
        return out

    def densities_from_log_ne(
        self, log_ne: float, log_a_cum: Sequence[float]
    ) -> list[float]:
        """Build [n_Ar, n_Ar+, n_Ar2+, n_Ar3+, n_Ar4+, n_e] from log(ne)."""
        ne = exp_clamped(log_ne)
        r1 = exp_clamped(log_a_cum[0] - log_ne)
        r2 = exp_clamped(log_a_cum[1] - 2.0 * log_ne)
        r3 = exp_clamped(log_a_cum[2] - 3.0 * log_ne)
        r4 = exp_clamped(log_a_cum[3] - 4.0 * log_ne)

        d_charge = r1 + 2.0 * r2 + 3.0 * r3 + 4.0 * r4
        d_charge = max(d_charge, MIN_POSITIVE_N)

        n0 = ne / d_charge
        n1 = n0 * r1
        n2 = n0 * r2
        n3 = n0 * r3
        n4 = n0 * r4
        return [n0, n1, n2, n3, n4, ne]

    def chemical_potentials_j(self, x_ln_n: Sequence[float], temperature_k: float) -> list[float]:
        mu: list[float] = []
        for idx, sp in enumerate(SPECIES):
            ln_n = x_ln_n[idx]
            lam = thermal_de_broglie_wavelength(self.masses[sp], temperature_k)
            q = self.q_int(sp, temperature_k)
            mu_i = K_B * temperature_k * (ln_n + 3.0 * math.log(lam) - math.log(q)) + self.eion_j[sp]
            mu.append(mu_i)
        return mu

    def residuals(self, x_ln_n: Sequence[float], temperature_k: float, pressure_pa: float) -> list[float]:
        n = [exp_clamped(x) for x in x_ln_n]
        mu = self.chemical_potentials_j(x_ln_n=x_ln_n, temperature_k=temperature_k)
        kt = K_B * temperature_k

        res = [0.0] * 6
        # Saha equilibrium relations in chemical-potential form
        # Scale by kT for dimensionless residuals.
        res[0] = (mu[0] - mu[1] - mu[5]) / kt  # Ar <-> Ar+ + e
        res[1] = (mu[1] - mu[2] - mu[5]) / kt  # Ar+ <-> Ar2+ + e
        res[2] = (mu[2] - mu[3] - mu[5]) / kt  # Ar2+ <-> Ar3+ + e
        res[3] = (mu[3] - mu[4] - mu[5]) / kt  # Ar3+ <-> Ar4+ + e

        # Quasi-neutrality in log-ratio form:
        # ne = n+  <=>  log(ne) - log(n+) = 0
        # This avoids Jacobian collapse when charged species are extremely small.
        n_charge = n[1] + 2.0 * n[2] + 3.0 * n[3] + 4.0 * n[4]
        res[4] = math.log(max(n[5], MIN_POSITIVE_N)) - math.log(max(n_charge, MIN_POSITIVE_N))

        # Equation of state with Debye-Huckel pressure correction
        p_ideal = sum(ni * kt for ni in n)
        p_dh = debye_huckel_pressure_correction(n=n, temperature_k=temperature_k)
        # Scale by pressure.
        res[5] = ((p_ideal - p_dh) - pressure_pa) / pressure_pa
        return res

    def analytic_jacobian(
        self, x_ln_n: Sequence[float], temperature_k: float, pressure_pa: float
    ) -> list[list[float]]:
        """Analytic Jacobian for the 6x6 residual system."""
        n = [exp_clamped(x) for x in x_ln_n]
        jac = [[0.0 for _ in range(6)] for _ in range(6)]

        # Saha (chemical potential) equations:
        # r0 = mu_Ar - mu_Ar+ - mu_e
        # r1 = mu_Ar+ - mu_Ar2+ - mu_e
        # r2 = mu_Ar2+ - mu_Ar3+ - mu_e
        # r3 = mu_Ar3+ - mu_Ar4+ - mu_e
        # Residuals are divided by kT -> coefficients are +/-1.
        for r in range(4):
            jac[r][r] += 1.0
            jac[r][r + 1] -= 1.0
            jac[r][5] -= 1.0

        # Quasi-neutrality in log-ratio form:
        # r4 = log(ne) - log(n_charge), n_charge = n1 + 2n2 + 3n3 + 4n4
        kt = K_B * temperature_k
        p_ref = pressure_pa
        n_charge = n[1] + 2.0 * n[2] + 3.0 * n[3] + 4.0 * n[4]
        n_charge_safe = max(n_charge, MIN_POSITIVE_N)
        jac[4][0] = 0.0
        jac[4][1] = -n[1] / n_charge_safe
        jac[4][2] = -2.0 * n[2] / n_charge_safe
        jac[4][3] = -3.0 * n[3] / n_charge_safe
        jac[4][4] = -4.0 * n[4] / n_charge_safe
        jac[4][5] = 1.0

        # EOS with Debye-Huckel: r5 = [sum(n_i kT) - deltaP_DH - P] / P_ref
        sigma = 0.0
        for ni, zi_abs in zip(n, ABS_CHARGES):
            sigma += (zi_abs * zi_abs) * ni
        delta_p_dh = debye_huckel_pressure_correction(n=n, temperature_k=temperature_k)

        for j in range(6):
            d_pideal = kt * n[j]
            if sigma > 0.0 and ABS_CHARGES[j] > 0:
                # deltaP = A * sigma^(3/2) -> d(deltaP)/dx_j
                #         = 1.5 * deltaP/sigma * (z_j^2 n_j)
                d_sigma = (ABS_CHARGES[j] ** 2) * n[j]
                d_dh = 1.5 * delta_p_dh / sigma * d_sigma
            else:
                d_dh = 0.0
            jac[5][j] = (d_pideal - d_dh) / p_ref

        return jac

    def solve_state_newton(
        self, temperature_k: float, pressure_pa: float, x0_ln_n: Sequence[float], tol: float, max_iter: int
    ) -> tuple[list[float], int, float]:
        x = list(x0_ln_n)
        for it in range(1, max_iter + 1):
            f = self.residuals(x_ln_n=x, temperature_k=temperature_k, pressure_pa=pressure_pa)
            fnorm = l2_norm(f)
            if fnorm < tol:
                return x, it, fnorm

            jac = self.analytic_jacobian(
                x_ln_n=x, temperature_k=temperature_k, pressure_pa=pressure_pa
            )
            dx = solve_linear_regularized(jac, [-v for v in f])

            # Damped update for robustness
            accepted = False
            alpha = 1.0
            for _ in range(20):
                x_trial = [
                    min(MAX_LN_N, max(MIN_LN_N, xi + alpha * dxi))
                    for xi, dxi in zip(x, dx)
                ]
                f_trial = self.residuals(
                    x_ln_n=x_trial, temperature_k=temperature_k, pressure_pa=pressure_pa
                )
                if l2_norm(f_trial) < fnorm:
                    x = x_trial
                    accepted = True
                    break
                alpha *= 0.5
            if not accepted:
                raise RuntimeError(
                    f"Newton line-search failed at T={temperature_k} K, P={pressure_pa} Pa."
                )

        final_res = l2_norm(self.residuals(x_ln_n=x, temperature_k=temperature_k, pressure_pa=pressure_pa))
        raise RuntimeError(
            f"Newton did not converge in {max_iter} iterations at T={temperature_k} K, "
            f"P={pressure_pa} Pa. Final residual norm={final_res:.3e}"
        )

    def solve_state_reduced_saha(
        self, temperature_k: float, pressure_pa: float, x0_ln_n: Sequence[float], tol: float, max_iter: int
    ) -> tuple[list[float], int, float]:
        """Fallback 1D solver: enforce Saha + charge exactly, solve EOS for ne."""
        kt = K_B * temperature_k
        n_ref = max(pressure_pa / kt, MIN_POSITIVE_N)

        log_k = self.saha_log_constants(temperature_k=temperature_k)
        log_a_cum = [
            log_k[0],
            log_k[0] + log_k[1],
            log_k[0] + log_k[1] + log_k[2],
            log_k[0] + log_k[1] + log_k[2] + log_k[3],
        ]

        def eos_res(log_ne: float) -> float:
            n = self.densities_from_log_ne(log_ne=log_ne, log_a_cum=log_a_cum)
            p_ideal = sum(ni * kt for ni in n)
            p_dh = debye_huckel_pressure_correction(n=n, temperature_k=temperature_k)
            return ((p_ideal - p_dh) - pressure_pa) / pressure_pa

        lo = max(MIN_LN_N, math.log(n_ref) - 120.0)
        hi = min(MAX_LN_N, math.log(n_ref) + 10.0)
        if hi <= lo:
            hi = min(MAX_LN_N, lo + 50.0)

        n_scan = 240
        xs = [lo + (hi - lo) * i / n_scan for i in range(n_scan + 1)]
        fs = [eos_res(x) for x in xs]

        bracket = None
        best_idx = 0
        best_abs = float("inf")
        for i, fi in enumerate(fs):
            afi = abs(fi)
            if math.isfinite(afi) and afi < best_abs:
                best_abs = afi
                best_idx = i
        for i in range(n_scan):
            f0 = fs[i]
            f1 = fs[i + 1]
            if not (math.isfinite(f0) and math.isfinite(f1)):
                continue
            if f0 == 0.0:
                bracket = (xs[i], xs[i], f0, f0)
                break
            if f0 * f1 < 0.0:
                bracket = (xs[i], xs[i + 1], f0, f1)
                break

        iters = 0
        if bracket is not None:
            a, b, fa, fb = bracket
            x = a if abs(fa) < abs(fb) else b
            for it in range(1, max_iter * 8 + 1):
                iters = it
                x = 0.5 * (a + b)
                fx = eos_res(x)
                if abs(fx) < tol:
                    break
                if fa * fx <= 0.0:
                    b, fb = x, fx
                else:
                    a, fa = x, fx
                if abs(b - a) < 1e-10:
                    break
        else:
            x = xs[best_idx]
            for it in range(1, max_iter * 8 + 1):
                iters = it
                fx = eos_res(x)
                if abs(fx) < tol:
                    break
                h = 1e-4
                xp = min(hi, x + h)
                xm = max(lo, x - h)
                fp = eos_res(xp)
                fm = eos_res(xm)
                dfdx = (fp - fm) / max(xp - xm, 1e-12)
                if not math.isfinite(dfdx) or abs(dfdx) < 1e-12:
                    step = -0.5 * fx
                else:
                    step = -fx / dfdx
                x_new = min(hi, max(lo, x + step))
                if abs(x_new - x) < 1e-10:
                    x = x_new
                    break
                x = x_new

        n = self.densities_from_log_ne(log_ne=x, log_a_cum=log_a_cum)
        if (not all(math.isfinite(v) for v in n)) or n[0] > 1e6 * n_ref:
            # Under-resolved ultra-low-ionization fallback in finite precision.
            n = [n_ref, 0.0, 0.0, 0.0, 0.0, 0.0]
        x_out = [ln_clamped_from_n(ni) for ni in n]
        fnorm = l2_norm(self.residuals(x_ln_n=x_out, temperature_k=temperature_k, pressure_pa=pressure_pa))
        return x_out, max(iters, 1), fnorm

    def solve_state(
        self, temperature_k: float, pressure_pa: float, x0_ln_n: Sequence[float], tol: float, max_iter: int
    ) -> tuple[list[float], int, float, str]:
        if scipy_root is not None:
            try:
                sol = scipy_root(
                    fun=lambda x: self.residuals(x_ln_n=x, temperature_k=temperature_k, pressure_pa=pressure_pa),
                    x0=list(x0_ln_n),
                    method="hybr",
                    options={"maxfev": max_iter * 20},
                )
                if sol.success:
                    fnorm = l2_norm(
                        self.residuals(
                            x_ln_n=sol.x.tolist(), temperature_k=temperature_k, pressure_pa=pressure_pa
                        )
                    )
                    return sol.x.tolist(), int(sol.nfev), fnorm, "scipy.root(hybr)"
            except Exception:
                pass

        try:
            x, it, fnorm = self.solve_state_newton(
                temperature_k=temperature_k,
                pressure_pa=pressure_pa,
                x0_ln_n=x0_ln_n,
                tol=tol,
                max_iter=max_iter,
            )
            return x, it, fnorm, "custom_newton"
        except RuntimeError:
            x, it, fnorm = self.solve_state_reduced_saha(
                temperature_k=temperature_k,
                pressure_pa=pressure_pa,
                x0_ln_n=x0_ln_n,
                tol=tol,
                max_iter=max_iter,
            )
            return x, it, fnorm, "reduced_saha_scalar"


def initial_guess_ln_n(temperature_k: float, pressure_pa: float) -> list[float]:
    n_total = pressure_pa / (K_B * temperature_k)
    n0 = max(0.999999 * n_total, 1e-300)
    n1 = max(1e-8 * n_total, 1e-300)
    n2 = max(1e-14 * n_total, 1e-300)
    n3 = max(1e-20 * n_total, 1e-300)
    n4 = max(1e-26 * n_total, 1e-300)
    ne = max(n1 + 2.0 * n2 + 3.0 * n3 + 4.0 * n4, 1e-300)
    return [math.log(v) for v in [n0, n1, n2, n3, n4, ne]]


def write_results_csv(out_csv: Path, rows: list[dict[str, float | str]]) -> None:
    fieldnames = [
        "T_K",
        "P_atm",
        "solver_used",
        "iterations",
        "residual_norm",
        "n_Ar",
        "n_Ar_p",
        "n_Ar2_p",
        "n_Ar3_p",
        "n_Ar4_p",
        "n_e",
        "charge_residual",
        "eos_residual_Pa",
        "deltaP_DH_Pa",
    ]
    with out_csv.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    out_dir = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    solver = ArgonLTEEquilibriumSolver(partition_json_path=args.partition_json.resolve())
    temperatures = build_temperature_grid(
        t_min=args.t_min,
        t_max=args.t_max,
        t_step=args.t_step,
        include_t_max=not args.exclude_t_max,
    )

    rows: list[dict[str, float | str]] = []
    run_log: dict[str, object] = {
        "solver_input": {
            "partition_json": str(args.partition_json.resolve()),
            "temperatures_K": {
                "t_min": temperatures[0],
                "t_max": temperatures[-1],
                "n_points": len(temperatures),
                "step": args.t_step,
            },
            "pressures_atm": args.pressures_atm,
            "tol": args.tol,
            "max_iter": args.max_iter,
        },
        "constants": {
            "k_B_J_per_K": K_B,
            "h_J_s": H_PLANCK,
            "epsilon0_F_per_m": EPS0,
            "e_C": E_CHARGE,
            "atm_to_pa": ATM_TO_PA,
        },
        "masses_kg": {
            "Ar": MASS_AR,
            "Ar+": MASS_AR - MASS_E,
            "Ar2+": MASS_AR - 2.0 * MASS_E,
            "Ar3+": MASS_AR - 3.0 * MASS_E,
            "Ar4+": MASS_AR - 4.0 * MASS_E,
            "e-": MASS_E,
        },
        "cumulative_ionization_eV_by_species": {
            sp: solver.eion_j[sp] / EV_TO_J for sp in SPECIES
        },
    }

    for p_atm in args.pressures_atm:
        p_pa = p_atm * ATM_TO_PA
        x_guess = initial_guess_ln_n(temperature_k=temperatures[0], pressure_pa=p_pa)
        for temperature_k in temperatures:
            x_sol, iters, fnorm, solver_used = solver.solve_state(
                temperature_k=temperature_k,
                pressure_pa=p_pa,
                x0_ln_n=x_guess,
                tol=args.tol,
                max_iter=args.max_iter,
            )
            x_guess = x_sol[:]  # warm-start next temperature
            n = [exp_clamped(v) for v in x_sol]
            charge_res = n[5] - (n[1] + 2.0 * n[2] + 3.0 * n[3] + 4.0 * n[4])
            p_ideal = sum(ni * K_B * temperature_k for ni in n)
            delta_p_dh = debye_huckel_pressure_correction(n=n, temperature_k=temperature_k)
            eos_res = (p_ideal - delta_p_dh) - p_pa

            rows.append(
                {
                    "T_K": temperature_k,
                    "P_atm": p_atm,
                    "solver_used": solver_used,
                    "iterations": iters,
                    "residual_norm": fnorm,
                    "n_Ar": n[0],
                    "n_Ar_p": n[1],
                    "n_Ar2_p": n[2],
                    "n_Ar3_p": n[3],
                    "n_Ar4_p": n[4],
                    "n_e": n[5],
                    "charge_residual": charge_res,
                    "eos_residual_Pa": eos_res,
                    "deltaP_DH_Pa": delta_p_dh,
                }
            )

    out_csv = out_dir / "argon_lte_equilibrium_0p1_1_4atm.csv"
    out_json = out_dir / "argon_lte_equilibrium_0p1_1_4atm.json"
    write_results_csv(out_csv=out_csv, rows=rows)
    out_json.write_text(
        json.dumps({"metadata": run_log, "rows": rows}, ensure_ascii=False),
        encoding="utf-8",
    )

    print(f"[OK] Wrote: {out_csv}")
    print(f"[OK] Wrote: {out_json}")
    print(f"[INFO] states solved: {len(rows)}")
    print(f"[INFO] pressures (atm): {args.pressures_atm}")
    print(f"[INFO] temperatures: {temperatures[0]} .. {temperatures[-1]} K ({len(temperatures)} points)")


if __name__ == "__main__":
    main()
