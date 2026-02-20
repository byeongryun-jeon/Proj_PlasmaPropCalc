#!/usr/bin/env python3
"""Build collision-integral transport input tables for Ar LTE plasma (Phase 4).

Implemented interactions:
    1) Ar-Ar      (neutral-neutral): LJ-style surrogate using Aziz1990 sigma, epsilon/k
    2) Ar-Ar+     (ion-neutral): resonant charge exchange + Langevin-like elastic
    3) e-Ar       (electron-neutral): Maxwellian average of Qm(E) from Milloy + Frost
    4) charged-charged: screened Coulomb using Mason1967 reduced tables + Debye length

Outputs:
    - non-charged collision-integral table (temperature-dependent)
    - charged collision-integral table (temperature + pressure-dependent)
    - merged table
    - metadata JSON with formulas/assumptions/source files
    - Mutation++ pair XML snippet for Ar/Ar+/e- related explicit pairs
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Callable


# --- Constants ---
K_B = 1.380649e-23  # J/K
K_B_EV = 8.617333262145e-5  # eV/K
EPS0 = 8.8541878128e-12  # F/m
E_CHARGE = 1.602176634e-19  # C
EV_TO_J = E_CHARGE
ANG2_TO_M2 = 1.0e-20
M2_TO_ANG2 = 1.0e20
PI = math.pi

# Ar static polarizability volume (for ion-induced dipole model)
ALPHA_AR_M3 = 1.6411e-30

# Charged species used in argon LTE phase-2/3 pipeline
CHARGED_SPECIES = ["e-", "Ar+", "Ar2+", "Ar3+", "Ar4+"]
CHARGES = {"e-": -1, "Ar+": 1, "Ar2+": 2, "Ar3+": 3, "Ar4+": 4}
EQ_CHARGE_COLS = {
    "e-": "n_e",
    "Ar+": "n_Ar_p",
    "Ar2+": "n_Ar2_p",
    "Ar3+": "n_Ar3_p",
    "Ar4+": "n_Ar4_p",
}

# Debye-Huckel reduced-temperature table (Mutation++ CoulombIntegrals.cpp).
# Values are (T*)^2 * Q^(1,4) and (T*)^2 * Q^(1,5) for attractive/repulsive
# screened Coulomb interactions. These are required by the 3rd-order electron
# Lee matrix in Devoto-style transport.
COULOMB_TSTAR_GRID = [
    0.1,
    0.2,
    0.3,
    0.4,
    0.6,
    0.8,
    1.0,
    2.0,
    3.0,
    4.0,
    6.0,
    8.0,
    10.0,
    20.0,
    30.0,
    40.0,
    60.0,
    80.0,
    100.0,
    200.0,
    300.0,
    400.0,
    600.0,
    800.0,
    1.0e3,
    1.0e4,
]

COULOMB_Q14_ATT_TSTAR2 = [
    0.0285,
    0.0460,
    0.0578,
    0.0669,
    0.0806,
    0.0910,
    0.0993,
    0.1269,
    0.1440,
    0.1566,
    0.1748,
    0.1880,
    0.1983,
    0.2307,
    0.2497,
    0.2634,
    0.2828,
    0.2969,
    0.3081,
    0.3433,
    0.3634,
    0.3778,
    0.3981,
    0.4125,
    0.4236,
    0.5388,
]

COULOMB_Q14_REP_TSTAR2 = [
    0.0110,
    0.0221,
    0.0316,
    0.0399,
    0.0536,
    0.0648,
    0.0742,
    0.1065,
    0.1268,
    0.1416,
    0.1627,
    0.1777,
    0.1894,
    0.2254,
    0.2462,
    0.2610,
    0.2817,
    0.2964,
    0.3078,
    0.3429,
    0.3633,
    0.3778,
    0.3981,
    0.4125,
    0.4236,
    0.5388,
]

COULOMB_Q15_ATT_TSTAR2 = [
    0.0227,
    0.0353,
    0.0437,
    0.0500,
    0.0596,
    0.0668,
    0.0725,
    0.0915,
    0.1032,
    0.1118,
    0.1241,
    0.1330,
    0.1400,
    0.1616,
    0.1744,
    0.1835,
    0.1967,
    0.2062,
    0.2137,
    0.2373,
    0.2506,
    0.2602,
    0.2737,
    0.2833,
    0.2908,
    0.3675,
]

COULOMB_Q15_REP_TSTAR2 = [
    0.0093,
    0.0181,
    0.0255,
    0.0317,
    0.0419,
    0.0500,
    0.0567,
    0.0792,
    0.0930,
    0.1030,
    0.1172,
    0.1272,
    0.1349,
    0.1589,
    0.1727,
    0.1825,
    0.1963,
    0.2061,
    0.2136,
    0.2370,
    0.2506,
    0.2602,
    0.2737,
    0.2833,
    0.2908,
    0.3675,
]

# Ar+-Ar resonant charge exchange constants (Murphy & Tam 2014)
AR_ARP_CX_STATES = [
    {"state": "2Sigma", "weight": 1.0 / 3.0, "A": 8.921, "B": 0.3960},
    {"state": "2Pi", "weight": 2.0 / 3.0, "A": 6.189, "B": 0.2934},
]


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Build Ar collision-integral transport inputs (Phase 4)."
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=root / "data" / "raw" / "murphy" / "collision_integrals" / "argon",
        help="Directory containing raw collision-integral source CSV files.",
    )
    parser.add_argument(
        "--equilibrium-csv",
        type=Path,
        default=root / "data" / "processed" / "equilibrium" / "argon_lte_equilibrium_0p1_1_4atm.csv",
        help="Phase-2 equilibrium CSV for electron density (charged-charged Debye length).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "data" / "processed" / "transport",
        help="Output directory for transport input tables.",
    )
    parser.add_argument("--t-min", type=float, default=300.0)
    parser.add_argument("--t-max", type=float, default=30000.0)
    parser.add_argument("--t-step", type=float, default=200.0)
    parser.add_argument(
        "--exclude-t-max",
        action="store_true",
        help="Do not force include t-max when off-grid.",
    )
    parser.add_argument(
        "--pressures-atm",
        type=float,
        nargs="+",
        default=[0.1, 1.0, 4.0],
        help="Pressure cases (atm) to include for charged-charged table.",
    )
    parser.add_argument(
        "--skip-charged",
        action="store_true",
        help="Skip charged-charged table generation.",
    )
    parser.add_argument(
        "--min-ne-m3",
        type=float,
        default=1.0e14,
        help="Minimum electron density used in Debye-length evaluation.",
    )
    parser.add_argument(
        "--max-debye-length-m",
        type=float,
        default=1.0e-3,
        help="Maximum Debye length cap for screened-Coulomb stability.",
    )
    parser.add_argument(
        "--debye-length-scale",
        type=float,
        default=1.0,
        help="Multiplier applied to Debye length before Coulomb reduced-temperature cutoff.",
    )
    parser.add_argument(
        "--coulomb-tstar-max",
        type=float,
        default=1.0e4,
        help="Upper cutoff of reduced temperature T* used in screened-Coulomb tables.",
    )
    parser.add_argument(
        "--coulomb-lnlambda-scale",
        type=float,
        default=1.0,
        help=(
            "Global scale factor for charged-charged Coulomb collision integrals "
            "(acts like a Coulomb-log convention adjustment)."
        ),
    )
    parser.add_argument(
        "--e-ar-high-order-blend",
        type=float,
        default=0.0,
        help=(
            "Blend factor for e-Ar Q14/Q15: 0 keeps Murphy closure (Q14=Q15=Q11), "
            "1 uses direct high-order moment integration."
        ),
    )
    parser.add_argument(
        "--debye-screening-density",
        choices=["electron_only", "all_charged"],
        default="electron_only",
        help=(
            "Screening density model for Debye length: "
            "electron_only uses n_e, all_charged uses sum(z_i^2 n_i) over charged species."
        ),
    )
    parser.add_argument(
        "--skip-mutationpp-snippet",
        action="store_true",
        help="Skip generating Mutation++ pair XML snippet.",
    )
    return parser.parse_args()


def build_temperature_grid(t_min: float, t_max: float, t_step: float, include_t_max: bool) -> list[float]:
    if t_step <= 0.0:
        raise ValueError("t-step must be > 0")
    if t_min <= 0.0:
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


def linear_interp(x_grid: list[float], y_grid: list[float], x: float) -> float:
    if x <= x_grid[0]:
        return y_grid[0]
    if x >= x_grid[-1]:
        return y_grid[-1]
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
    return y0 + w * (y1 - y0)


def log_interp(x_grid: list[float], y_grid: list[float], x: float) -> float:
    x_clipped = min(max(x, x_grid[0]), x_grid[-1])
    lx_grid = [math.log10(v) for v in x_grid]
    return linear_interp(lx_grid, y_grid, math.log10(x_clipped))


def debye_huckel_lambda_tstar(
    temperature_k: float,
    electron_density_m3: float,
    zabs: int,
    max_debye_length_m: float,
    debye_length_scale: float,
    coulomb_tstar_max: float,
) -> tuple[float, float]:
    """Return Debye length and reduced temperature for screened Coulomb CI."""
    t = max(temperature_k, 1.0e-12)
    ne = max(electron_density_m3, 1.0e-16)

    # LTE diffusion notes use electron-density-only Debye screening.
    lambda_d = math.sqrt(EPS0 * K_B * t / (ne * E_CHARGE * E_CHARGE))
    lambda_d *= max(debye_length_scale, 1.0e-12)

    z_eff = max(float(zabs), 1.0)
    b = z_eff * E_CHARGE * E_CHARGE / (8.0 * PI * EPS0 * K_B * t)

    # Keep T* inside the tabulated range, consistent with Debye-Huckel tables.
    tstar_max = max(coulomb_tstar_max, COULOMB_TSTAR_GRID[0])
    lambda_cap_tstar = 2.0 * tstar_max * b
    lambda_d = min(lambda_d, max_debye_length_m, lambda_cap_tstar)

    t_star = max(0.5 * lambda_d / max(b, 1.0e-300), COULOMB_TSTAR_GRID[0])
    return lambda_d, t_star


def trapz(x: list[float], y: list[float]) -> float:
    if len(x) != len(y):
        raise ValueError("x/y length mismatch")
    total = 0.0
    for i in range(len(x) - 1):
        total += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    return total


def thermal_average_q(
    temperature_k: float,
    q_of_energy_a2: Callable[[float], float],
    moment_power: int,
    n_grid: int = 800,
) -> float:
    """Return Maxwellian average of Q(E) in Angstrom^2 over energy E in eV."""
    if temperature_k <= 0.0:
        raise ValueError("temperature must be > 0")

    kT_eV = K_B_EV * temperature_k
    e_min = 1.0e-4
    e_max = max(80.0, 40.0 * kT_eV)
    if e_max <= e_min:
        e_max = e_min * 10.0

    log_min = math.log(e_min)
    log_max = math.log(e_max)
    energies = [
        math.exp(log_min + (log_max - log_min) * i / (n_grid - 1))
        for i in range(n_grid)
    ]

    numerand: list[float] = []
    denomand: list[float] = []
    for e_ev in energies:
        w = (e_ev**moment_power) * math.exp(-e_ev / kT_eV)
        q = max(q_of_energy_a2(e_ev), 0.0)
        numerand.append(q * w)
        denomand.append(w)

    den = trapz(energies, denomand)
    if den <= 0.0:
        return 0.0
    num = trapz(energies, numerand)
    return num / den


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        return [dict(row) for row in reader]


def load_single_constant(path: Path, param_name: str) -> float:
    for row in read_csv_rows(path):
        if row.get("param") == param_name:
            return float(row["value"])
    raise ValueError(f"Parameter {param_name!r} not found in {path}")


def load_e_ar_qm_data(milloy_path: Path, frost_path: Path) -> tuple[list[float], list[float], list[str]]:
    milloy = read_csv_rows(milloy_path)
    frost = read_csv_rows(frost_path)

    data: dict[float, tuple[float, str]] = {}
    for row in milloy:
        e = float(row["E_eV"])
        q = float(row["Qm_A2"])
        data[e] = (q, row.get("source", "Milloy1977"))

    # Use Frost extension for E > max(Milloy energy)
    e_max_milloy = max(data.keys())
    for row in frost:
        e = float(row["E_eV"])
        if e <= e_max_milloy:
            continue
        q = float(row["Qm_A2"])
        data[e] = (q, row.get("source", "Frost1964"))

    energies = sorted(data.keys())
    qvals = [data[e][0] for e in energies]
    srcs = [data[e][1] for e in energies]
    return energies, qvals, srcs


def build_piecewise_q_interpolator(
    energies_eV: list[float], qvals_a2: list[float]
) -> Callable[[float], float]:
    if len(energies_eV) != len(qvals_a2):
        raise ValueError("Energy/Q arrays mismatch")
    if len(energies_eV) < 2:
        raise ValueError("Need at least 2 data points")

    e_grid = energies_eV
    q_grid = qvals_a2

    def q_func(e_ev: float) -> float:
        e = max(e_ev, 1.0e-8)
        if e <= e_grid[0]:
            return q_grid[0]
        if e >= e_grid[-1]:
            # Smooth high-energy tail for numerical stability.
            return q_grid[-1] * math.sqrt(e_grid[-1] / e)
        return linear_interp(e_grid, q_grid, e)

    return q_func


def lj_omega11_star(t_star: float) -> float:
    t = max(t_star, 1.0e-10)
    return (
        1.06036 / (t**0.15610)
        + 0.19300 * math.exp(-0.47635 * t)
        + 1.03587 * math.exp(-1.52996 * t)
        + 1.76474 * math.exp(-3.89411 * t)
    )


def lj_omega22_star(t_star: float) -> float:
    t = max(t_star, 1.0e-10)
    return (
        1.16145 / (t**0.14874)
        + 0.52487 * math.exp(-0.77320 * t)
        + 2.16178 * math.exp(-2.43787 * t)
    )


def q_ar_ar_m2(temperature_k: float, epsilon_over_k: float, sigma_angstrom: float) -> tuple[float, float]:
    t_star = max(temperature_k / epsilon_over_k, 1.0e-12)
    omega11_star = lj_omega11_star(t_star)
    omega22_star = lj_omega22_star(t_star)
    sigma_m = sigma_angstrom * 1.0e-10
    area = PI * sigma_m * sigma_m
    return omega11_star * area, omega22_star * area


def q_ar_arp_charge_exchange_a2(energy_ev: float) -> float:
    e = max(energy_ev, 1.0e-10)
    ln_e = math.log(e)
    q_ex = 0.0
    for state in AR_ARP_CX_STATES:
        q_ex += state["weight"] * ((state["A"] - state["B"] * ln_e) ** 2)
    return 2.0 * q_ex


def q_ar_arp_elastic_langevin_a2(energy_ev: float, charge_state: int = 1) -> float:
    """Capture-model elastic cross-section for V(r) = -C4/r^4 in Angstrom^2."""
    e_j = max(energy_ev * EV_TO_J, 1.0e-30)
    z = float(abs(charge_state))
    # Polarizability in this project is stored as a volume [m^3].
    # Convert to SI polarizability alpha_SI = 4*pi*eps0*alpha_vol before
    # evaluating C4 in SI. Without this conversion, sigma is overestimated
    # by ~O(1e5) in the Ar-Ar+ elastic capture model.
    alpha_si = 4.0 * PI * EPS0 * ALPHA_AR_M3
    c4 = alpha_si * (z * E_CHARGE) ** 2 / (32.0 * PI * PI * EPS0 * EPS0)
    sigma_m2 = PI * math.sqrt(2.0 * c4 / e_j)
    return sigma_m2 * M2_TO_ANG2


def load_mason_tables(
    table_i_path: Path, table_ii_path: Path
) -> tuple[dict[str, list[float]], dict[str, list[float]]]:
    rows_i = read_csv_rows(table_i_path)
    rows_ii = read_csv_rows(table_ii_path)

    ti = {
        "Tstar": [float(r["Tstar"]) for r in rows_i],
        "O11_attr_t2": [float(r["Omega11_t2_attraction"]) for r in rows_i],
        "O11_rep_t2": [float(r["Omega11_t2_repulsion"]) for r in rows_i],
        "O22_attr_t2": [float(r["Omega22_t2_attraction"]) for r in rows_i],
        "O22_rep_t2": [float(r["Omega22_t2_repulsion"]) for r in rows_i],
    }
    tii = {
        "Tstar": [float(r["Tstar"]) for r in rows_ii],
        "A_attr": [float(r["Astar_attraction"]) for r in rows_ii],
        "A_rep": [float(r["Astar_repulsion"]) for r in rows_ii],
        "B_attr": [float(r["Bstar_attraction"]) for r in rows_ii],
        "B_rep": [float(r["Bstar_repulsion"]) for r in rows_ii],
        "C_attr": [float(r["Cstar_attraction"]) for r in rows_ii],
        "C_rep": [float(r["Cstar_repulsion"]) for r in rows_ii],
    }
    return ti, tii


def load_equilibrium_screening_density_by_pressure(
    path: Path,
    mode: str,
) -> dict[float, tuple[list[float], list[float]]]:
    rows = read_csv_rows(path)
    grouped: dict[float, list[tuple[float, float]]] = {}
    for r in rows:
        p = float(r["P_atm"])
        t = float(r["T_K"])
        ne = max(float(r["n_e"]), 0.0)

        if mode == "electron_only":
            screening_density = ne
        else:
            screening_density = 0.0
            for sp, col in EQ_CHARGE_COLS.items():
                ni = max(float(r[col]), 0.0)
                zi = abs(CHARGES[sp])
                screening_density += (zi * zi) * ni

        grouped.setdefault(p, []).append((t, screening_density))

    out: dict[float, tuple[list[float], list[float]]] = {}
    for p, arr in grouped.items():
        arr_sorted = sorted(arr, key=lambda x: x[0])
        out[p] = ([v[0] for v in arr_sorted], [v[1] for v in arr_sorted])
    return out


def nearest_pressure_key(available: list[float], target: float, tol: float = 1.0e-10) -> float:
    for p in available:
        if abs(p - target) <= tol:
            return p
    raise ValueError(f"Requested pressure {target} atm not found. Available: {available}")


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def format_xml_series(values: list[float], precision: int = 6) -> str:
    return " ".join(f"{v:.{precision}g}" for v in values)


def format_xml_temperature_nodes(values: list[float]) -> str:
    return " ".join(str(int(round(v))) for v in values)


def generate_mutationpp_snippet(
    out_path: Path,
    temp_nodes: list[float],
    non_charged_rows: list[dict[str, object]],
) -> None:
    def pair_series(pair: str, key: str) -> list[float]:
        pair_rows = [r for r in non_charged_rows if r["pair"] == pair]
        t_grid = [float(r["T_K"]) for r in pair_rows]
        y_grid = [float(r[key]) for r in pair_rows]
        return [linear_interp(t_grid, y_grid, t) for t in temp_nodes]

    q11_e_ar = pair_series("e-Ar", "Q11_A2")
    q22_e_ar = pair_series("e-Ar", "Q22_A2")
    q11_ar_ar = pair_series("Ar-Ar", "Q11_A2")
    q22_ar_ar = pair_series("Ar-Ar", "Q22_A2")
    q11_ar_arp = pair_series("Ar-Ar+", "Q11_A2")
    q22_ar_arp = pair_series("Ar-Ar+", "Q22_A2")

    text = (
        "<!-- Auto-generated by build_collision_integrals.py -->\n"
        "<!-- Approximate Phase-4 transport inputs for Ar LTE tutorial -->\n"
        "<collision-pairs>\n"
        "  <pair s1=\"e-\" s2=\"Ar\">\n"
        "    <Q11 type=\"table\" units=\"K,Å-Å\" ref=\"Milloy1977+Frost1964\">\n"
        f"      {format_xml_temperature_nodes(temp_nodes)},\n"
        f"      {format_xml_series(q11_e_ar, 5)}\n"
        "    </Q11>\n"
        "    <Q22 type=\"table\" units=\"K,Å-Å\" ref=\"Murphy2014 assumption Omega22=Omega11\">\n"
        f"      {format_xml_temperature_nodes(temp_nodes)},\n"
        f"      {format_xml_series(q22_e_ar, 5)}\n"
        "    </Q22>\n"
        "    <Bst type=\"constant\" value=\"1.0\" ref=\"placeholder\"/>\n"
        "    <Cst type=\"constant\" value=\"1.0\" ref=\"placeholder\"/>\n"
        "  </pair>\n\n"
        "  <pair s1=\"Ar\" s2=\"Ar\">\n"
        "    <Q11 type=\"table\" units=\"K,Å-Å\" ref=\"Aziz1990+LJ surrogate\">\n"
        f"      {format_xml_temperature_nodes(temp_nodes)},\n"
        f"      {format_xml_series(q11_ar_ar, 5)}\n"
        "    </Q11>\n"
        "    <Q22 type=\"table\" units=\"K,Å-Å\" ref=\"Aziz1990+LJ surrogate\">\n"
        f"      {format_xml_temperature_nodes(temp_nodes)},\n"
        f"      {format_xml_series(q22_ar_ar, 5)}\n"
        "    </Q22>\n"
        "    <Bst type=\"constant\" value=\"1.15\" ref=\"Wright2005 default neutral-neutral\"/>\n"
        "    <Cst type=\"constant\" value=\"0.92\" ref=\"Wright2005 default neutral-neutral\"/>\n"
        "  </pair>\n\n"
        "  <pair s1=\"Ar\" s2=\"Ar+\">\n"
        "    <Q11 type=\"table\" units=\"K,Å-Å\" ref=\"MurphyTam2014+Aubreton1986 model\">\n"
        f"      {format_xml_temperature_nodes(temp_nodes)},\n"
        f"      {format_xml_series(q11_ar_arp, 5)}\n"
        "    </Q11>\n"
        "    <Q22 type=\"table\" units=\"K,Å-Å\" ref=\"MurphyTam2014+Aubreton1986 model\">\n"
        f"      {format_xml_temperature_nodes(temp_nodes)},\n"
        f"      {format_xml_series(q22_ar_arp, 5)}\n"
        "    </Q22>\n"
        "    <Bst type=\"constant\" value=\"1.20\" ref=\"Wright2005 default ion-neutral\"/>\n"
        "    <Cst type=\"constant\" value=\"0.85\" ref=\"Wright2005 default ion-neutral\"/>\n"
        "  </pair>\n"
        "</collision-pairs>\n"
    )
    out_path.write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.raw_dir.resolve()

    temperatures = build_temperature_grid(
        t_min=args.t_min,
        t_max=args.t_max,
        t_step=args.t_step,
        include_t_max=not args.exclude_t_max,
    )

    # --- Load raw sources ---
    aziz_path = raw_dir / "aziz1990_hfdtcs2_constants.csv"
    milloy_path = raw_dir / "milloy1977_e_ar_momentum_transfer_table1.csv"
    frost_path = raw_dir / "frost1964_e_ar_high_energy_digitized.csv"
    mason_i_path = raw_dir / "mason1967_tableI_t2omega.csv"
    mason_ii_path = raw_dir / "mason1967_tableII_abc.csv"
    aubreton_path = raw_dir / "aubreton1986_ar_arplus_elastic_potentials.csv"

    epsilon_over_k = load_single_constant(aziz_path, "epsilon_over_k")
    sigma_angstrom = load_single_constant(aziz_path, "sigma")
    e_grid, q_grid, _ = load_e_ar_qm_data(milloy_path, frost_path)
    q_e_ar = build_piecewise_q_interpolator(e_grid, q_grid)
    mason_i, mason_ii = load_mason_tables(mason_i_path, mason_ii_path)

    # --- Non-charged categories ---
    non_charged_rows: list[dict[str, object]] = []
    for t in temperatures:
        # Ar-Ar (neutral-neutral)
        q11_m2, q22_m2 = q_ar_ar_m2(t, epsilon_over_k, sigma_angstrom)
        ast = q22_m2 / max(q11_m2, 1.0e-300)
        non_charged_rows.append(
            {
                "T_K": t,
                "pair": "Ar-Ar",
                "category": "neutral-neutral",
                "Q11_m2": q11_m2,
                "Q22_m2": q22_m2,
                "Q11_A2": q11_m2 * M2_TO_ANG2,
                "Q22_A2": q22_m2 * M2_TO_ANG2,
                "Ast": ast,
                "Bst": 1.15,
                "Cst": 0.92,
                "source": "Aziz1990_TableI_II + LJ surrogate",
                "notes": "HFDTCS2 full scattering integral not yet implemented",
            }
        )

        # Ar-Ar+ (ion-neutral): combine inelastic and elastic with sqrt(sum squares)
        q_in_11 = thermal_average_q(t, q_ar_arp_charge_exchange_a2, moment_power=1)
        q_in_22 = thermal_average_q(t, q_ar_arp_charge_exchange_a2, moment_power=2)
        q_el_11 = thermal_average_q(
            t, lambda e: q_ar_arp_elastic_langevin_a2(e, charge_state=1), moment_power=1
        )
        q_el_22 = thermal_average_q(
            t, lambda e: q_ar_arp_elastic_langevin_a2(e, charge_state=1), moment_power=2
        )
        q11_ar_arp_a2 = math.sqrt(q_in_11 * q_in_11 + q_el_11 * q_el_11)
        q22_ar_arp_a2 = math.sqrt(q_in_22 * q_in_22 + q_el_22 * q_el_22)
        q11_m2 = q11_ar_arp_a2 * ANG2_TO_M2
        q22_m2 = q22_ar_arp_a2 * ANG2_TO_M2
        ast = q22_m2 / max(q11_m2, 1.0e-300)
        non_charged_rows.append(
            {
                "T_K": t,
                "pair": "Ar-Ar+",
                "category": "ion-neutral",
                "Q11_m2": q11_m2,
                "Q22_m2": q22_m2,
                "Q11_A2": q11_ar_arp_a2,
                "Q22_A2": q22_ar_arp_a2,
                "Ast": ast,
                "Bst": 1.20,
                "Cst": 0.85,
                "source": "MurphyTam2014 charge-exchange + Langevin elastic",
                "notes": "Uses Aubreton1986 context; elastic represented with capture model",
            }
        )

        # Ar-Ar(z+) for z=2,3,4:
        # keep explicit rows in the collision table using the same scaling
        # convention previously applied as runtime fallback in Phase 5.
        for z in (2, 3, 4):
            scale_z = math.sqrt(float(z))
            q11_z_a2 = q11_ar_arp_a2 * scale_z
            q22_z_a2 = q22_ar_arp_a2 * scale_z
            q11_z_m2 = q11_z_a2 * ANG2_TO_M2
            q22_z_m2 = q22_z_a2 * ANG2_TO_M2
            ast_z = q22_z_m2 / max(q11_z_m2, 1.0e-300)
            non_charged_rows.append(
                {
                    "T_K": t,
                    "pair": f"Ar-Ar{z}+",
                    "category": "ion-neutral-derived",
                    "Q11_m2": q11_z_m2,
                    "Q22_m2": q22_z_m2,
                    "Q11_A2": q11_z_a2,
                    "Q22_A2": q22_z_a2,
                    "Ast": ast_z,
                    "Bst": 1.20,
                    "Cst": 0.85,
                    "source": "Scaled from Ar-Ar+",
                    "notes": "Q11/Q22 scaled by sqrt(z) from Ar-Ar+ for higher ionization states",
                }
            )

        # e-Ar (electron-neutral)
        q11_a2 = thermal_average_q(t, q_e_ar, moment_power=1)
        q22_a2 = q11_a2  # Murphy-style approximation for e-neutral
        # Devoto electron Lee matrix needs high-order moments (Q14/Q15).
        # Use regularized moment closure: blend direct high-order moments with
        # Murphy closure to keep low-temperature conductivity stable.
        q14_a2_moment = thermal_average_q(t, q_e_ar, moment_power=4)
        q15_a2_moment = thermal_average_q(t, q_e_ar, moment_power=5)
        blend = min(max(args.e_ar_high_order_blend, 0.0), 1.0)
        q14_a2 = (1.0 - blend) * q11_a2 + blend * q14_a2_moment
        q15_a2 = (1.0 - blend) * q11_a2 + blend * q15_a2_moment
        q11_m2 = q11_a2 * ANG2_TO_M2
        q22_m2 = q22_a2 * ANG2_TO_M2
        q14_m2 = q14_a2 * ANG2_TO_M2
        q15_m2 = q15_a2 * ANG2_TO_M2
        non_charged_rows.append(
            {
                "T_K": t,
                "pair": "e-Ar",
                "category": "electron-neutral",
                "Q11_m2": q11_m2,
                "Q22_m2": q22_m2,
                "Q14_m2": q14_m2,
                "Q15_m2": q15_m2,
                "Q11_A2": q11_a2,
                "Q22_A2": q22_a2,
                "Q14_A2": q14_a2,
                "Q15_A2": q15_a2,
                "Ast": 1.0,
                "Bst": 1.0,
                "Cst": 1.0,
                "source": "Milloy1977 + Frost1964",
                "notes": "Q14/Q15 from regularized direct Maxwellian moments of Qm(E); Omega22 assumed equal to Omega11",
            }
        )

    # --- Charged-charged category ---
    charged_rows: list[dict[str, object]] = []
    if not args.skip_charged:
        eq_by_p = load_equilibrium_screening_density_by_pressure(
            args.equilibrium_csv.resolve(),
            mode=args.debye_screening_density,
        )
        available_p = sorted(eq_by_p.keys())

        tstar_grid_i = mason_i["Tstar"]
        tstar_grid_ii = mason_ii["Tstar"]

        for p_target in args.pressures_atm:
            p_key = nearest_pressure_key(available_p, p_target)
            t_eq, screening_eq = eq_by_p[p_key]
            for t in temperatures:
                n_screen_raw = max(linear_interp(t_eq, screening_eq, t), 1.0e-300)
                n_screen = max(n_screen_raw, args.min_ne_m3)
                for i in range(len(CHARGED_SPECIES)):
                    for j in range(i, len(CHARGED_SPECIES)):
                        s1 = CHARGED_SPECIES[i]
                        s2 = CHARGED_SPECIES[j]
                        z1 = CHARGES[s1]
                        z2 = CHARGES[s2]
                        zabs = abs(z1 * z2)
                        if zabs == 0:
                            continue

                        lambda_d, t_star = debye_huckel_lambda_tstar(
                            temperature_k=t,
                            electron_density_m3=n_screen,
                            zabs=zabs,
                            max_debye_length_m=args.max_debye_length_m,
                            debye_length_scale=args.debye_length_scale,
                            coulomb_tstar_max=args.coulomb_tstar_max,
                        )

                        interaction = "attraction" if z1 * z2 < 0 else "repulsion"

                        if interaction == "attraction":
                            o11_t2 = log_interp(tstar_grid_i, mason_i["O11_attr_t2"], t_star)
                            o22_t2 = log_interp(tstar_grid_i, mason_i["O22_attr_t2"], t_star)
                            a_star = log_interp(tstar_grid_ii, mason_ii["A_attr"], t_star)
                            b_star = log_interp(tstar_grid_ii, mason_ii["B_attr"], t_star)
                            c_star = log_interp(tstar_grid_ii, mason_ii["C_attr"], t_star)
                            q14_t2 = log_interp(COULOMB_TSTAR_GRID, COULOMB_Q14_ATT_TSTAR2, t_star)
                            q15_t2 = log_interp(COULOMB_TSTAR_GRID, COULOMB_Q15_ATT_TSTAR2, t_star)
                        else:
                            o11_t2 = log_interp(tstar_grid_i, mason_i["O11_rep_t2"], t_star)
                            o22_t2 = log_interp(tstar_grid_i, mason_i["O22_rep_t2"], t_star)
                            a_star = log_interp(tstar_grid_ii, mason_ii["A_rep"], t_star)
                            b_star = log_interp(tstar_grid_ii, mason_ii["B_rep"], t_star)
                            c_star = log_interp(tstar_grid_ii, mason_ii["C_rep"], t_star)
                            q14_t2 = log_interp(COULOMB_TSTAR_GRID, COULOMB_Q14_REP_TSTAR2, t_star)
                            q15_t2 = log_interp(COULOMB_TSTAR_GRID, COULOMB_Q15_REP_TSTAR2, t_star)

                        scale = PI * lambda_d * lambda_d / max(t_star * t_star, 1.0e-300)
                        coulomb_scale = max(args.coulomb_lnlambda_scale, 1.0e-12)
                        q11_m2 = coulomb_scale * o11_t2 * scale
                        q22_m2 = coulomb_scale * o22_t2 * scale
                        q14_m2 = coulomb_scale * q14_t2 * scale
                        q15_m2 = coulomb_scale * q15_t2 * scale

                        charged_rows.append(
                            {
                                "T_K": t,
                                "P_atm": p_key,
                                "pair": f"{s1}__{s2}",
                                "interaction": interaction,
                                "lambda_D_m": lambda_d,
                                "Tstar": t_star,
                                "Q11_m2": q11_m2,
                                "Q22_m2": q22_m2,
                                "Q14_m2": q14_m2,
                                "Q15_m2": q15_m2,
                                "Q11_A2": q11_m2 * M2_TO_ANG2,
                                "Q22_A2": q22_m2 * M2_TO_ANG2,
                                "Q14_A2": q14_m2 * M2_TO_ANG2,
                                "Q15_A2": q15_m2 * M2_TO_ANG2,
                                "Ast": a_star,
                                "Bst": b_star,
                                "Cst": c_star,
                                "source": "Mason1967 table + Debye-Huckel screening",
                            }
                        )

    # --- Write outputs ---
    non_charged_csv = output_dir / "argon_collision_integrals_non_charged.csv"
    non_charged_fields = [
        "T_K",
        "pair",
        "category",
        "Q11_m2",
        "Q22_m2",
        "Q14_m2",
        "Q15_m2",
        "Q11_A2",
        "Q22_A2",
        "Q14_A2",
        "Q15_A2",
        "Ast",
        "Bst",
        "Cst",
        "source",
        "notes",
    ]
    write_csv(non_charged_csv, non_charged_fields, non_charged_rows)

    charged_csv = output_dir / "argon_collision_integrals_charged_by_pressure.csv"
    charged_fields = [
        "T_K",
        "P_atm",
        "pair",
        "interaction",
        "lambda_D_m",
        "Tstar",
        "Q11_m2",
        "Q22_m2",
        "Q14_m2",
        "Q15_m2",
        "Q11_A2",
        "Q22_A2",
        "Q14_A2",
        "Q15_A2",
        "Ast",
        "Bst",
        "Cst",
        "source",
    ]
    if charged_rows:
        write_csv(charged_csv, charged_fields, charged_rows)

    all_csv = output_dir / "argon_collision_integrals_all.csv"
    all_fields = [
        "T_K",
        "P_atm",
        "pair",
        "category",
        "interaction",
        "lambda_D_m",
        "Tstar",
        "Q11_m2",
        "Q22_m2",
        "Q14_m2",
        "Q15_m2",
        "Q11_A2",
        "Q22_A2",
        "Q14_A2",
        "Q15_A2",
        "Ast",
        "Bst",
        "Cst",
        "source",
        "notes",
    ]
    all_rows: list[dict[str, object]] = []
    for r in non_charged_rows:
        all_rows.append(
            {
                "T_K": r["T_K"],
                "P_atm": "all",
                "pair": r["pair"],
                "category": r["category"],
                "interaction": "",
                "lambda_D_m": "",
                "Tstar": "",
                "Q11_m2": r["Q11_m2"],
                "Q22_m2": r["Q22_m2"],
                "Q14_m2": r.get("Q14_m2", ""),
                "Q15_m2": r.get("Q15_m2", ""),
                "Q11_A2": r["Q11_A2"],
                "Q22_A2": r["Q22_A2"],
                "Q14_A2": r.get("Q14_A2", ""),
                "Q15_A2": r.get("Q15_A2", ""),
                "Ast": r["Ast"],
                "Bst": r["Bst"],
                "Cst": r["Cst"],
                "source": r["source"],
                "notes": r["notes"],
            }
        )
    for r in charged_rows:
        all_rows.append(
            {
                "T_K": r["T_K"],
                "P_atm": r["P_atm"],
                "pair": r["pair"],
                "category": "charged-charged",
                "interaction": r["interaction"],
                "lambda_D_m": r["lambda_D_m"],
                "Tstar": r["Tstar"],
                "Q11_m2": r["Q11_m2"],
                "Q22_m2": r["Q22_m2"],
                "Q14_m2": r["Q14_m2"],
                "Q15_m2": r["Q15_m2"],
                "Q11_A2": r["Q11_A2"],
                "Q22_A2": r["Q22_A2"],
                "Q14_A2": r["Q14_A2"],
                "Q15_A2": r["Q15_A2"],
                "Ast": r["Ast"],
                "Bst": r["Bst"],
                "Cst": r["Cst"],
                "source": r["source"],
                "notes": "",
            }
        )
    all_rows.sort(key=lambda r: (str(r["pair"]), str(r["P_atm"]), float(r["T_K"])))
    write_csv(all_csv, all_fields, all_rows)

    snippet_path = None
    if not args.skip_mutationpp_snippet:
        snippet_path = output_dir / "mutationpp_argon_collision_pairs.xml"
        temp_nodes = [1000.0, 2000.0, 4000.0, 5000.0, 6000.0, 8000.0, 10000.0, 15000.0, 20000.0]
        generate_mutationpp_snippet(
            out_path=snippet_path, temp_nodes=temp_nodes, non_charged_rows=non_charged_rows
        )

    metadata = {
        "inputs": {
            "raw_dir": str(raw_dir),
            "equilibrium_csv": str(args.equilibrium_csv.resolve()),
            "source_files": {
                "aziz1990": str(aziz_path),
                "aubreton1986": str(aubreton_path),
                "milloy1977": str(milloy_path),
                "frost1964_digitized": str(frost_path),
                "mason1967_tableI": str(mason_i_path),
                "mason1967_tableII": str(mason_ii_path),
            },
        },
        "temperature_grid": {
            "t_min": min(temperatures),
            "t_max": max(temperatures),
            "t_step": args.t_step,
            "n_points": len(temperatures),
        },
        "pressures_atm_for_charged": args.pressures_atm if not args.skip_charged else [],
        "charged_model_limits": {
            "min_ne_m3": args.min_ne_m3,
            "max_debye_length_m": args.max_debye_length_m,
            "debye_length_scale": args.debye_length_scale,
            "coulomb_tstar_max": args.coulomb_tstar_max,
            "coulomb_lnlambda_scale": args.coulomb_lnlambda_scale,
            "e_ar_high_order_blend": args.e_ar_high_order_blend,
        },
        "debye_screening_density_model": args.debye_screening_density,
        "models": {
            "Ar-Ar": "LJ surrogate with Aziz epsilon/k and sigma",
            "Ar-Ar+": "Qex model (MurphyTam2014) + Langevin elastic capture model",
            "Ar-Ar2+/Ar-Ar3+/Ar-Ar4+": "Derived from Ar-Ar+ via sqrt(z) scaling on Q11/Q22",
            "e-Ar": "Maxwellian average of Qm(E), Milloy table + Frost high-energy extension",
            "charged-charged": "Mason1967 reduced tables + Debye-Huckel screening for pairs among e-, Ar+, Ar2+, Ar3+, Ar4+",
        },
        "outputs": {
            "non_charged_csv": str(non_charged_csv),
            "charged_csv": str(charged_csv if charged_rows else ""),
            "all_csv": str(all_csv),
            "mutationpp_snippet_xml": str(snippet_path) if snippet_path is not None else "",
        },
        "notes": [
            "Ar-Ar full HFDTCS2 classical scattering integration is planned; current version uses a documented LJ-style surrogate.",
            "Frost1964 high-energy e-Ar points are figure-digitized approximations.",
            "Charged-charged table uses electron density from phase-2 equilibrium CSV and is pressure-dependent.",
            "Debye length is capped and ne is floored for low-ionization stability in screened-Coulomb mode.",
            "Coulomb tuning knobs (debye_length_scale, coulomb_tstar_max, coulomb_lnlambda_scale) can be calibrated against reference lookup tables.",
        ],
    }
    metadata_path = output_dir / "argon_collision_integrals_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"[OK] Wrote: {non_charged_csv}")
    if charged_rows:
        print(f"[OK] Wrote: {charged_csv}")
    print(f"[OK] Wrote: {all_csv}")
    print(f"[OK] Wrote: {metadata_path}")
    if snippet_path is not None:
        print(f"[OK] Wrote: {snippet_path}")
    print(f"[INFO] Temperature points: {len(temperatures)}")
    if not args.skip_charged:
        print(f"[INFO] Charged rows: {len(charged_rows)}")


if __name__ == "__main__":
    main()
