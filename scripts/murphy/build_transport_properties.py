#!/usr/bin/env python3
"""Build Phase-5 transport coefficients (mu, kappa, sigma) for Ar LTE plasma.

Inputs:
    - Phase 2 equilibrium composition CSV (n_i(T,P))
    - Phase 3 thermo CSV (rho, h, cp)
    - Phase 4 collision integrals CSV (long or wide schema)
    - Phase 1 partition-function JSON (species frozen enthalpy/cp)

Outputs:
    - combined CSV and per-pressure CSV transport tables
    - metadata JSON (mapping, assumptions, fallback usage, peak summaries)
    - validation plots (mu, kappa total, sigma; plus kappa component splits)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


# --- Constants ---
K_B = 1.380649e-23  # J/K
K_B_EV = 8.617333262145e-5  # eV/K
E_CHARGE = 1.602176634e-19  # C
PI = math.pi
EPS0 = 8.8541878128e-12  # F/m
SQRT2 = math.sqrt(2.0)

MASS_AR = 6.6335209e-26
MASS_E = 9.1093837e-31

SPECIES = ["Ar", "Ar+", "Ar2+", "Ar3+", "Ar4+", "e-"]
HEAVY = ["Ar", "Ar+", "Ar2+", "Ar3+", "Ar4+"]
SPECIES_ORDER = {sp: i for i, sp in enumerate(SPECIES)}

MASS = {
    "Ar": MASS_AR,
    "Ar+": MASS_AR - MASS_E,
    "Ar2+": MASS_AR - 2.0 * MASS_E,
    "Ar3+": MASS_AR - 3.0 * MASS_E,
    "Ar4+": MASS_AR - 4.0 * MASS_E,
    "e-": MASS_E,
}

CHARGE = {
    "Ar": 0,
    "Ar+": 1,
    "Ar2+": 2,
    "Ar3+": 3,
    "Ar4+": 4,
    "e-": -1,
}

# Debye-Huckel reduced-temperature table (Mutation++ CoulombIntegrals.cpp).
# Stored values correspond to:
#   (T*)^2 * Q22_rep, (T*)^2 * Q24_rep, and E*_rep = Q23_rep / Q22_rep
# for repulsive Coulomb interactions (e-e).
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

COULOMB_Q22_REP_TSTAR2 = [
    0.0304,
    0.0697,
    0.1086,
    0.1459,
    0.2144,
    0.2757,
    0.3310,
    0.5460,
    0.6999,
    0.8197,
    1.0006,
    1.1385,
    1.2435,
    1.5892,
    1.7959,
    1.9438,
    2.1531,
    2.3019,
    2.4171,
    2.7713,
    2.9747,
    3.1177,
    3.3185,
    3.4610,
    3.5719,
    4.7211,
]

COULOMB_Q24_REP_TSTAR2 = [
    0.0208,
    0.0445,
    0.0661,
    0.0856,
    0.1193,
    0.1476,
    0.1719,
    0.2587,
    0.3154,
    0.3574,
    0.4182,
    0.4620,
    0.4962,
    0.6027,
    0.6649,
    0.7089,
    0.7707,
    0.8145,
    0.8483,
    0.9534,
    1.0149,
    1.0583,
    1.1193,
    1.1624,
    1.1959,
    1.5413,
]

COULOMB_EST_REP = [
    0.8146,
    0.7836,
    0.7634,
    0.7483,
    0.7265,
    0.7108,
    0.6988,
    0.6628,
    0.6436,
    0.6311,
    0.6152,
    0.6052,
    0.5981,
    0.5797,
    0.5714,
    0.5663,
    0.5600,
    0.5561,
    0.5533,
    0.5456,
    0.5419,
    0.5398,
    0.5373,
    0.5358,
    0.5348,
    0.5265,
]

EQ_COL_TO_SPECIES = {
    "n_Ar": "Ar",
    "n_Ar_p": "Ar+",
    "n_Ar2_p": "Ar2+",
    "n_Ar3_p": "Ar3+",
    "n_Ar4_p": "Ar4+",
    "n_e": "e-",
}


@dataclass
class CollisionRecord:
    t_k: float
    p_atm: float | None
    s1: str
    s2: str
    q11_m2: float
    q22_m2: float
    q14_m2: float | None
    q15_m2: float | None
    ast: float
    bst: float
    cst: float
    category: str
    source: str


@dataclass
class CollisionValue:
    q11_m2: float
    q22_m2: float
    q14_m2: float | None
    q15_m2: float | None
    ast: float
    bst: float
    cst: float
    category: str
    source: str
    fallback: str


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Build Phase-5 transport coefficients (mu, kappa, sigma)."
    )
    parser.add_argument(
        "--equilibrium-csv",
        type=Path,
        default=root / "data" / "processed" / "equilibrium" / "argon_lte_equilibrium_0p1_1_4atm.csv",
        help="Phase-2 equilibrium composition CSV.",
    )
    parser.add_argument(
        "--thermo-csv",
        type=Path,
        default=root / "data" / "processed" / "thermo" / "argon_thermo_0p1_1_4atm.csv",
        help="Phase-3 thermo CSV.",
    )
    parser.add_argument(
        "--collision-csv",
        type=Path,
        default=root / "data" / "processed" / "transport" / "argon_collision_integrals_all.csv",
        help="Phase-4 collision integral CSV (long or wide format).",
    )
    parser.add_argument(
        "--partition-json",
        type=Path,
        default=root / "data" / "processed" / "partition_functions" / "partition_functions_Ar_system.json",
        help="Phase-1 partition-function JSON.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "data" / "processed" / "transport_properties",
        help="Output directory for Phase-5 transport tables.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip gnuplot figure generation.",
    )
    parser.add_argument(
        "--matf-reference-csv",
        type=Path,
        default=root / "data" / "processed" / "thermo" / "argon_matf_reference_0p1_1_4atm.csv",
        help="Optional MATF reference CSV for dashed overlay on summary plots.",
    )
    parser.add_argument(
        "--reaction-scale",
        type=float,
        default=1.0,
        help="Multiplier for k_reac term (default: 1.0).",
    )
    parser.add_argument(
        "--pr-heavy",
        type=float,
        default=0.67,
        help="Heavy translational Prandtl-like constant.",
    )
    parser.add_argument(
        "--pr-int",
        type=float,
        default=1.0,
        help="Internal conductivity Prandtl-like constant.",
    )
    parser.add_argument(
        "--k-reac-model",
        choices=["legacy_composite", "butler_mass", "butler_mole", "butler_mole_cp"],
        default="legacy_composite",
        help=(
            "Reaction-conductivity closure: "
            "legacy_composite=(mass-gradient + cp amplification), "
            "butler_mass=(mass-gradient only), "
            "butler_mole=(mole-gradient only), "
            "butler_mole_cp=(mole-gradient + cp amplification)."
        ),
    )
    parser.add_argument(
        "--murphy-strict",
        action="store_true",
        help=(
            "Apply stricter Murphy-style transport closure defaults. "
            "Currently sets k_reac model to Butler mole-gradient form."
        ),
    )
    return parser.parse_args()


def apply_murphy_profile(args: argparse.Namespace) -> None:
    if not args.murphy_strict:
        return
    args.k_reac_model = "butler_mole"
    args.reaction_scale = 1.0


def canonical_pair(s1: str, s2: str) -> tuple[str, str]:
    if SPECIES_ORDER[s1] <= SPECIES_ORDER[s2]:
        return s1, s2
    return s2, s1


def derive_bst_cst(ast: float) -> tuple[float, float]:
    a = max(ast, 1.0e-8)
    bst = 1.0 + 0.35 * (a - 1.0)
    cst = 0.75 + 0.25 / a
    return bst, cst


def numeric_derivative(x: list[float], y: list[float]) -> list[float]:
    if len(x) != len(y):
        raise ValueError("x/y length mismatch")
    n = len(x)
    if n < 2:
        raise ValueError("Need at least two points")
    out = [0.0] * n

    dx0 = x[1] - x[0]
    if abs(dx0) < 1e-30:
        dx0 = 1e-30
    out[0] = (y[1] - y[0]) / dx0

    dxn = x[-1] - x[-2]
    if abs(dxn) < 1e-30:
        dxn = 1e-30
    out[-1] = (y[-1] - y[-2]) / dxn

    for i in range(1, n - 1):
        dx = x[i + 1] - x[i - 1]
        if abs(dx) < 1e-30:
            dx = 1e-30
        out[i] = (y[i + 1] - y[i - 1]) / dx
    return out


def interp1(x: list[float], y: list[float], xq: float) -> float:
    if xq <= x[0]:
        return y[0]
    if xq >= x[-1]:
        return y[-1]
    lo, hi = 0, len(x) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x[mid] <= xq:
            lo = mid
        else:
            hi = mid
    x0, x1 = x[lo], x[hi]
    y0, y1 = y[lo], y[hi]
    w = (xq - x0) / (x1 - x0)
    return y0 + w * (y1 - y0)


def parse_pair_from_string(pair: str) -> tuple[str, str] | None:
    pair = pair.strip()
    for s1 in SPECIES:
        for s2 in SPECIES:
            candidates = [f"{s1}__{s2}", f"{s1}-{s2}"]
            # Phase-4 CSV uses shorthand like "e-Ar" (not "e--Ar").
            if s1 == "e-":
                candidates.append(f"e-{s2}")
            if s2 == "e-":
                candidates.append(f"{s1}-e")
            if pair in candidates:
                return s1, s2
    return None


def detect_collision_format(headers: list[str]) -> str:
    h = set(headers)
    if "pair" in h and ("Q11_m2" in h or "Q11_A2" in h):
        return "long"
    pair_pat = re.compile(r".+_(?:Q)?[12]_[12]$")
    if any(pair_pat.match(col) for col in headers):
        return "wide"
    raise ValueError(
        "Unable to detect collision CSV schema. Expected long format with "
        "pair/Q11/Q22 columns or wide columns like Ar-Ar_1_1."
    )


def col_as_float(row: dict[str, str], keys: list[str], default: float | None = None) -> float:
    for key in keys:
        if key in row and row[key] not in ("", None):
            return float(row[key])
    if default is None:
        raise KeyError(f"None of keys found: {keys}")
    return default


def load_collision_long(path: Path, headers: list[str]) -> tuple[list[CollisionRecord], dict[str, object]]:
    q11_col = "Q11_m2" if "Q11_m2" in headers else "Q11_A2"
    q22_col = "Q22_m2" if "Q22_m2" in headers else "Q22_A2"
    q14_col = "Q14_m2" if "Q14_m2" in headers else ("Q14_A2" if "Q14_A2" in headers else "")
    q15_col = "Q15_m2" if "Q15_m2" in headers else ("Q15_A2" if "Q15_A2" in headers else "")
    area_scale = 1.0 if q11_col.endswith("_m2") else 1.0e-20

    records: list[CollisionRecord] = []
    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            if not row:
                continue
            pair_raw = row.get("pair", "")
            pair = parse_pair_from_string(pair_raw)
            if pair is None:
                continue

            t_k = float(row["T_K"])
            p_raw = row.get("P_atm", "")
            p_atm: float | None
            if p_raw in ("", "all", "ALL", "Any", "any"):
                p_atm = None
            else:
                p_atm = float(p_raw)

            q11 = float(row[q11_col]) * area_scale
            q22 = float(row[q22_col]) * area_scale
            if q11 <= 0.0:
                continue
            if q22 <= 0.0:
                q22 = q11

            q14: float | None = None
            q15: float | None = None
            if q14_col and row.get(q14_col, "") not in ("", None):
                v = float(row[q14_col]) * (1.0 if q14_col.endswith("_m2") else 1.0e-20)
                q14 = v if v > 0.0 else None
            if q15_col and row.get(q15_col, "") not in ("", None):
                v = float(row[q15_col]) * (1.0 if q15_col.endswith("_m2") else 1.0e-20)
                q15 = v if v > 0.0 else None

            ast = col_as_float(row, ["Ast", "Astar", "A_st"], default=q22 / q11)
            bst_val = row.get("Bst", "")
            cst_val = row.get("Cst", "")
            if bst_val in ("", None) or cst_val in ("", None):
                bst, cst = derive_bst_cst(ast)
            else:
                bst = float(bst_val)
                cst = float(cst_val)

            records.append(
                CollisionRecord(
                    t_k=t_k,
                    p_atm=p_atm,
                    s1=pair[0],
                    s2=pair[1],
                    q11_m2=q11,
                    q22_m2=q22,
                    q14_m2=q14,
                    q15_m2=q15,
                    ast=ast,
                    bst=bst,
                    cst=cst,
                    category=row.get("category", ""),
                    source=row.get("source", ""),
                )
            )

    mapping = {
        "format": "long",
        "q11_column": q11_col,
        "q22_column": q22_col,
        "q14_column": q14_col if q14_col else "",
        "q15_column": q15_col if q15_col else "",
        "ast_column": "Ast" if "Ast" in headers else "derived(Q22/Q11)",
        "bst_column": "Bst" if "Bst" in headers else "derived",
        "cst_column": "Cst" if "Cst" in headers else "derived",
    }
    return records, mapping


def parse_wide_collision_column(col: str) -> tuple[str, int, int] | None:
    m = re.match(r"(?P<pair>.+)_(?:Q)?(?P<l>[12])_(?P<s>[12])$", col)
    if not m:
        return None
    return m.group("pair"), int(m.group("l")), int(m.group("s"))


def load_collision_wide(path: Path, headers: list[str]) -> tuple[list[CollisionRecord], dict[str, object]]:
    t_col = "T_K" if "T_K" in headers else "T"
    p_col = "P_atm" if "P_atm" in headers else None

    mapped_cols: dict[str, dict[tuple[int, int], str]] = defaultdict(dict)
    for col in headers:
        parsed = parse_wide_collision_column(col)
        if parsed is None:
            continue
        pair, l, s = parsed
        mapped_cols[pair][(l, s)] = col

    records: list[CollisionRecord] = []
    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            t_k = float(row[t_col])
            p_atm = float(row[p_col]) if p_col is not None and row[p_col] != "" else None
            for pair_str, cols in mapped_cols.items():
                pair = parse_pair_from_string(pair_str)
                if pair is None:
                    continue
                col11 = cols.get((1, 1))
                col22 = cols.get((2, 2))
                if col11 is None:
                    continue
                q11 = float(row[col11])
                q22 = float(row[col22]) if col22 is not None and row[col22] != "" else q11
                if q11 <= 0.0:
                    continue
                ast = q22 / q11
                bst, cst = derive_bst_cst(ast)
                records.append(
                    CollisionRecord(
                        t_k=t_k,
                        p_atm=p_atm,
                        s1=pair[0],
                        s2=pair[1],
                        q11_m2=q11,
                        q22_m2=q22,
                        q14_m2=None,
                        q15_m2=None,
                        ast=ast,
                        bst=bst,
                        cst=cst,
                        category="",
                        source="wide-mapped",
                    )
                )

    mapping = {
        "format": "wide",
        "mapped_pairs": {
            pair: {
                "11": cols.get((1, 1), ""),
                "22": cols.get((2, 2), ""),
            }
            for pair, cols in mapped_cols.items()
        },
    }
    return records, mapping


def load_collision_records(path: Path) -> tuple[list[CollisionRecord], dict[str, object]]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.reader(fp)
        headers = next(reader)

    fmt = detect_collision_format(headers)
    if fmt == "long":
        return load_collision_long(path, headers)
    return load_collision_wide(path, headers)


class CollisionDatabase:
    def __init__(self, records: list[CollisionRecord]) -> None:
        self._tables: dict[tuple[str, str, float | None], dict[str, list[float] | str]] = {}
        self.fallback_counts: dict[str, int] = defaultdict(int)

        grouped: dict[tuple[str, str, float | None], list[CollisionRecord]] = defaultdict(list)
        for rec in records:
            s1, s2 = canonical_pair(rec.s1, rec.s2)
            p_key = None if rec.p_atm is None else float(rec.p_atm)
            grouped[(s1, s2, p_key)].append(rec)

        for key, rows in grouped.items():
            rows_sorted = sorted(rows, key=lambda r: r.t_k)
            table: dict[str, list[float] | str] = {
                "T": [r.t_k for r in rows_sorted],
                "Q11": [r.q11_m2 for r in rows_sorted],
                "Q22": [r.q22_m2 for r in rows_sorted],
                "Ast": [r.ast for r in rows_sorted],
                "Bst": [r.bst for r in rows_sorted],
                "Cst": [r.cst for r in rows_sorted],
                "category": rows_sorted[0].category,
                "source": rows_sorted[0].source,
            }
            if all(r.q14_m2 is not None for r in rows_sorted):
                table["Q14"] = [float(r.q14_m2) for r in rows_sorted if r.q14_m2 is not None]
            if all(r.q15_m2 is not None for r in rows_sorted):
                table["Q15"] = [float(r.q15_m2) for r in rows_sorted if r.q15_m2 is not None]
            self._tables[key] = table

    def _interp(self, key: tuple[str, str, float | None], t_k: float) -> CollisionValue | None:
        table = self._tables.get(key)
        if table is None:
            return None
        t = table["T"]
        q11 = interp1(t, table["Q11"], t_k)
        q22 = interp1(t, table["Q22"], t_k)
        q14 = interp1(t, table["Q14"], t_k) if "Q14" in table else None
        q15 = interp1(t, table["Q15"], t_k) if "Q15" in table else None
        ast = interp1(t, table["Ast"], t_k)
        bst = interp1(t, table["Bst"], t_k)
        cst = interp1(t, table["Cst"], t_k)
        cat = str(table["category"])
        src = str(table["source"])
        return CollisionValue(q11, q22, q14, q15, ast, bst, cst, cat, src, fallback="")

    def _lookup_direct(self, s1: str, s2: str, t_k: float, p_atm: float) -> CollisionValue | None:
        a, b = canonical_pair(s1, s2)
        val = self._interp((a, b, p_atm), t_k)
        if val is not None:
            return val
        return self._interp((a, b, None), t_k)

    def get(self, s1: str, s2: str, t_k: float, p_atm: float, track: bool = True) -> CollisionValue:
        direct = self._lookup_direct(s1, s2, t_k, p_atm)
        if direct is not None:
            return direct

        pair_set = {s1, s2}
        if "Ar" in pair_set:
            ion = s2 if s1 == "Ar" else s1
            if ion in ("Ar2+", "Ar3+", "Ar4+"):
                base = self._lookup_direct("Ar", "Ar+", t_k, p_atm)
                if base is not None:
                    z = abs(CHARGE[ion])
                    scale = math.sqrt(float(z))
                    if track:
                        self.fallback_counts[f"Ar-{ion}:scaled_from_Ar-Ar+"] += 1
                    return CollisionValue(
                        q11_m2=base.q11_m2 * scale,
                        q22_m2=base.q22_m2 * scale,
                        q14_m2=None,
                        q15_m2=None,
                        ast=base.ast,
                        bst=base.bst,
                        cst=base.cst,
                        category="ion-neutral-fallback",
                        source=f"scaled:{base.source}",
                        fallback="scaled_from_Ar-Ar+",
                    )

        self_i = self._lookup_direct(s1, s1, t_k, p_atm)
        self_j = self._lookup_direct(s2, s2, t_k, p_atm)
        if self_i is not None and self_j is not None:
            if track:
                self.fallback_counts[f"{s1}-{s2}:geometric_self"] += 1
            q11 = math.sqrt(max(self_i.q11_m2, 1e-300) * max(self_j.q11_m2, 1e-300))
            q22 = math.sqrt(max(self_i.q22_m2, 1e-300) * max(self_j.q22_m2, 1e-300))
            ast = q22 / max(q11, 1e-300)
            bst = 0.5 * (self_i.bst + self_j.bst)
            cst = 0.5 * (self_i.cst + self_j.cst)
            return CollisionValue(
                q11_m2=q11,
                q22_m2=q22,
                q14_m2=None,
                q15_m2=None,
                ast=ast,
                bst=bst,
                cst=cst,
                category="geometric-fallback",
                source=f"geometric({self_i.source},{self_j.source})",
                fallback="geometric_self",
            )

        if track:
            self.fallback_counts[f"{s1}-{s2}:constant"] += 1
        return CollisionValue(
            q11_m2=1.0e-19,
            q22_m2=1.2e-19,
            q14_m2=None,
            q15_m2=None,
            ast=1.2,
            bst=1.1,
            cst=0.9,
            category="constant-fallback",
            source="constant_fallback",
            fallback="constant",
        )


def load_equilibrium(path: Path) -> tuple[dict[float, list[dict[str, float]]], list[float], list[float]]:
    by_pressure: dict[float, list[dict[str, float]]] = defaultdict(list)
    temps: set[float] = set()

    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            p = float(row["P_atm"])
            t = float(row["T_K"])
            temps.add(t)
            state = {"T_K": t, "P_atm": p}
            for col, sp in EQ_COL_TO_SPECIES.items():
                state[sp] = float(row[col])
            by_pressure[p].append(state)

    pressures = sorted(by_pressure.keys())
    for p in pressures:
        by_pressure[p].sort(key=lambda r: r["T_K"])
    return by_pressure, sorted(temps), pressures


def load_thermo(path: Path) -> dict[tuple[float, float], dict[str, float]]:
    out: dict[tuple[float, float], dict[str, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            t = float(row["T_K"])
            p = float(row["P_atm"])
            out[(p, t)] = {
                "rho": float(row["rho_kg_m3"]),
                "h": float(row["h_J_kg"]),
                "cp": float(row["cp_J_kgK"]),
            }
    return out


def load_species_frozen_thermo(
    partition_json: Path,
) -> tuple[list[float], dict[str, list[float]], dict[str, list[float]], dict[str, float]]:
    payload = json.loads(partition_json.read_text(encoding="utf-8"))
    t_grid = [float(v) for v in payload["temperature_K"]]
    q = payload["Q"]
    eion = {sp: float(v) for sp, v in payload["cumulative_ionization_eV_by_species"].items()}

    h_mass: dict[str, list[float]] = {}
    cp_mass: dict[str, list[float]] = {}

    for sp in SPECIES:
        qv = [max(float(v), 1.0e-300) for v in q[sp]]
        lnq = [math.log(v) for v in qv]
        dlnq_dt = numeric_derivative(t_grid, lnq)

        h_ev_particle = [
            2.5 * K_B_EV * t + K_B_EV * (t**2) * dl + eion[sp]
            for t, dl in zip(t_grid, dlnq_dt)
        ]
        h_j_particle = [v * E_CHARGE for v in h_ev_particle]
        h_j_kg = [v / MASS[sp] for v in h_j_particle]
        cp_j_kgk = numeric_derivative(t_grid, h_j_kg)

        h_mass[sp] = h_j_kg
        cp_mass[sp] = cp_j_kgk

    return t_grid, h_mass, cp_mass, eion


def reduced_mass(m1: float, m2: float) -> float:
    return (m1 * m2) / max(m1 + m2, 1.0e-300)


def collision_rate_coeff(t_k: float, q_m2: float, m1: float, m2: float) -> float:
    mu = reduced_mass(m1, m2)
    v_th = math.sqrt(max(8.0 * K_B * t_k / (PI * mu), 0.0))
    return max(q_m2, 0.0) * v_th


def _gaussian_solve(a: list[list[float]], b: list[float]) -> list[float]:
    n = len(a)
    m = [row[:] for row in a]
    rhs = b[:]

    for k in range(n):
        piv = k
        piv_abs = abs(m[k][k])
        for i in range(k + 1, n):
            v = abs(m[i][k])
            if v > piv_abs:
                piv = i
                piv_abs = v
        if piv_abs < 1e-30:
            raise RuntimeError("Singular matrix")
        if piv != k:
            m[k], m[piv] = m[piv], m[k]
            rhs[k], rhs[piv] = rhs[piv], rhs[k]

        akk = m[k][k]
        for i in range(k + 1, n):
            factor = m[i][k] / akk
            if factor == 0.0:
                continue
            for j in range(k, n):
                m[i][j] -= factor * m[k][j]
            rhs[i] -= factor * rhs[k]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = rhs[i]
        for j in range(i + 1, n):
            s -= m[i][j] * x[j]
        aii = m[i][i]
        if abs(aii) < 1e-30:
            raise RuntimeError("Singular back substitution")
        x[i] = s / aii
    return x


def solve_linear_system(a: list[list[float]], b: list[float]) -> list[float]:
    n = len(a)
    diag_scale = max(abs(a[i][i]) for i in range(n))
    if diag_scale < 1.0:
        diag_scale = 1.0

    regs = [0.0, 1e-14, 1e-12, 1e-10, 1e-8]
    for reg in regs:
        aa = [row[:] for row in a]
        if reg > 0.0:
            for i in range(n):
                aa[i][i] += reg * diag_scale
        try:
            return _gaussian_solve(aa, b)
        except RuntimeError:
            continue
    # Last resort: diagonal fallback.
    x = [0.0] * n
    for i in range(n):
        den = a[i][i]
        if abs(den) < 1e-30:
            den = 1e-30
        x[i] = b[i] / den
    return x


def matrix_metrics(mat: list[list[float]]) -> dict[str, float]:
    n = len(mat)
    if n == 0:
        return {
            "diag_min": 0.0,
            "diag_max": 0.0,
            "diag_ratio": 0.0,
            "offdiag_to_diag": 0.0,
        }
    diag = [abs(mat[i][i]) for i in range(n)]
    diag_min = min(diag)
    diag_max = max(diag)
    diag_sum = sum(diag)
    offdiag_sum = 0.0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            offdiag_sum += abs(mat[i][j])
    return {
        "diag_min": diag_min,
        "diag_max": diag_max,
        "diag_ratio": diag_max / max(diag_min, 1.0e-300),
        "offdiag_to_diag": offdiag_sum / max(diag_sum, 1.0e-300),
    }


def normalize_weights(weights: dict[str, float]) -> dict[str, float]:
    total = sum(max(v, 0.0) for v in weights.values())
    if total <= 0.0:
        return {k: 0.0 for k in weights}
    return {k: max(v, 0.0) / total for k, v in weights.items()}


def dominant_pair(weight_share: dict[str, float]) -> tuple[str, float]:
    if not weight_share:
        return "", 0.0
    pair = max(weight_share, key=weight_share.get)
    return pair, weight_share[pair]


def debye_huckel_repulsive_ee_integrals(
    t_k: float, n_e: float
) -> tuple[float, float, float, float, float]:
    """Return (Q22ee, Q23ee, Q24ee, lambda_D, Tstar) for e-e repulsive Coulomb."""
    t = max(t_k, 1.0e-12)
    ne = max(n_e, 1.0e-16)

    b = E_CHARGE * E_CHARGE / (8.0 * PI * EPS0 * K_B * t)
    lambda_d = math.sqrt(EPS0 * K_B * t / (ne * E_CHARGE * E_CHARGE))
    lambda_capped = min(lambda_d, 2.0 * COULOMB_TSTAR_GRID[-1] * b)
    tstar = max(0.5 * lambda_capped / max(b, 1.0e-300), COULOMB_TSTAR_GRID[0])

    q22_tst2 = interp1(COULOMB_TSTAR_GRID, COULOMB_Q22_REP_TSTAR2, tstar)
    q24_tst2 = interp1(COULOMB_TSTAR_GRID, COULOMB_Q24_REP_TSTAR2, tstar)
    est = interp1(COULOMB_TSTAR_GRID, COULOMB_EST_REP, tstar)

    scale = PI * lambda_capped * lambda_capped / max(tstar * tstar, 1.0e-300)
    q22 = max(q22_tst2 * scale, 1.0e-300)
    q24 = max(q24_tst2 * scale, 1.0e-300)
    q23 = max(est * q22, 1.0e-300)
    return q22, q23, q24, lambda_capped, tstar


def weighted_collision_stars(
    species_list: list[str],
    n: dict[str, float],
    t_k: float,
    p_atm: float,
    db: CollisionDatabase,
) -> tuple[float, float, float]:
    wsum = 0.0
    a_acc = 0.0
    b_acc = 0.0
    c_acc = 0.0

    for i, si in enumerate(species_list):
        for sj in species_list[i:]:
            nij = n[si] * n[sj]
            if nij <= 0.0:
                continue
            val = db.get(si, sj, t_k, p_atm)
            k11 = collision_rate_coeff(t_k, val.q11_m2, MASS[si], MASS[sj])
            w = nij * k11
            wsum += w
            a_acc += w * val.ast
            b_acc += w * val.bst
            c_acc += w * val.cst

    if wsum <= 0.0:
        return 1.0, 1.0, 1.0
    return a_acc / wsum, b_acc / wsum, c_acc / wsum


def compute_transport_for_state(
    state: dict[str, float],
    thermo: dict[str, float],
    y_deriv: dict[str, float],
    x_deriv: dict[str, float],
    t_pf: list[float],
    h_mass_pf: dict[str, list[float]],
    cp_mass_pf: dict[str, list[float]],
    db: CollisionDatabase,
    pr_heavy: float,
    pr_int: float,
    reaction_scale: float,
    k_reac_model: str,
) -> tuple[dict[str, float], dict[str, object]]:
    t = state["T_K"]
    p = state["P_atm"]

    n = {sp: float(state[sp]) for sp in SPECIES}
    rho = max(float(thermo["rho"]), 1e-30)
    h_mix = float(thermo["h"])
    cp_mix = float(thermo["cp"])

    n_total = max(sum(n.values()), 1.0e-300)
    n_heavy = sum(n[sp] for sp in HEAVY)
    x = {sp: n[sp] / n_total for sp in SPECIES}
    y = {sp: (n[sp] * MASS[sp]) / rho for sp in SPECIES}

    nh = len(HEAVY)
    xh = [max(x[sp], 1.0e-16) for sp in HEAVY]

    # Precompute heavy-species collision data and factors used by Devoto matrices.
    etai: list[float] = [0.0] * nh
    for i, sp in enumerate(HEAVY):
        q22ii = max(db.get(sp, sp, t, p).q22_m2, 1.0e-300)
        etafac = 5.0 / 16.0 * math.sqrt(PI * K_B) * math.sqrt(MASS[sp])
        etai[i] = math.sqrt(t) * etafac / q22ii

    heavy_pair_weights_raw: dict[str, float] = {}
    for i, si in enumerate(HEAVY):
        for sj in HEAVY[i:]:
            val_ij = db.get(si, sj, t, p, track=False)
            k11_ij = collision_rate_coeff(t, val_ij.q11_m2, MASS[si], MASS[sj])
            heavy_pair_weights_raw[f"{si}-{sj}"] = n[si] * n[sj] * k11_ij
    heavy_pair_weights = normalize_weights(heavy_pair_weights_raw)

    pair_data: dict[tuple[int, int], tuple[float, float, float, float]] = {}
    for j in range(nh):
        for i in range(j + 1, nh):
            si = HEAVY[i]
            sj = HEAVY[j]
            val = db.get(si, sj, t, p)
            q11ij = max(val.q11_m2, 1.0e-300)
            dijfac = 3.0 / 16.0 * math.sqrt(
                2.0 * PI * K_B * (MASS[si] + MASS[sj]) / (MASS[si] * MASS[sj])
            )
            nDij = math.sqrt(t) * dijfac / q11ij
            pair_data[(i, j)] = (nDij, val.ast, val.bst, val.cst)

    # --- Viscosity mu (Devoto / CE 1st order heavy matrix) ---
    a_mu = [[0.0 for _ in range(nh)] for _ in range(nh)]
    for i in range(nh):
        a_mu[i][i] = xh[i] * xh[i] / max(etai[i], 1.0e-300)

    for j in range(nh):
        for i in range(j + 1, nh):
            nDij, ast, _, _ = pair_data[(i, j)]
            mi = MASS[HEAVY[i]]
            mj = MASS[HEAVY[j]]
            fac = xh[i] * xh[j] / max(nDij * (mi + mj), 1.0e-300)
            aij = fac * (1.2 * ast - 2.0)
            a_mu[i][j] = aij
            a_mu[j][i] = aij
            a_mu[i][i] += fac * (1.2 * mj / mi * ast + 2.0)
            a_mu[j][j] += fac * (1.2 * mi / mj * ast + 2.0)

    mu_matrix_stat = matrix_metrics(a_mu)
    alpha_mu = solve_linear_system(a_mu, xh)
    mu_heavy = max(sum(xh[i] * alpha_mu[i] for i in range(nh)), 1.0e-12)

    # --- Heavy translational conductivity k_H (Devoto / CE 2nd order matrix) ---
    a_h = [[0.0 for _ in range(nh)] for _ in range(nh)]
    fac0 = 4.0 / (15.0 * K_B)
    for i in range(nh):
        mi = MASS[HEAVY[i]]
        a_h[i][i] = fac0 * xh[i] * xh[i] * mi / max(etai[i], 1.0e-300)

    for j in range(nh):
        for i in range(j + 1, nh):
            nDij, ast, bst, _ = pair_data[(i, j)]
            mi = MASS[HEAVY[i]]
            mj = MASS[HEAVY[j]]
            mij_sum = mi + mj
            miij = mi / mij_sum
            mjij = mj / mij_sum
            fac = xh[i] * xh[j] / max(nDij * 25.0 * K_B, 1.0e-300)
            aij = fac * miij * mjij * (16.0 * ast + 12.0 * bst - 55.0)
            a_h[i][j] = aij
            a_h[j][i] = aij
            a_h[i][i] += fac * (
                miij * (30.0 * miij + 16.0 * mjij * ast) + mjij * mjij * (25.0 - 12.0 * bst)
            )
            a_h[j][j] += fac * (
                mjij * (30.0 * mjij + 16.0 * miij * ast) + miij * miij * (25.0 - 12.0 * bst)
            )

    kh_matrix_stat = matrix_metrics(a_h)
    alpha_h = solve_linear_system(a_h, xh)
    k_h = max(sum(xh[i] * alpha_h[i] for i in range(nh)), 0.0)

    # --- Electron subsystem: strict Devoto Lee matrix (order 3) for sigma and k_e ---
    q11_eh: list[float] = []
    q12_eh: list[float] = []
    q13_eh: list[float] = []
    q14_eh: list[float] = []
    q15_eh: list[float] = []
    wsum_e = 0.0
    a_e_acc = 0.0
    b_e_acc = 0.0
    c_e_acc = 0.0

    electron_pair_weights_raw: dict[str, float] = {}
    for sp in HEAVY:
        val = db.get("e-", sp, t, p)
        q11 = max(val.q11_m2, 1.0e-300)
        q12 = max(val.cst * q11, 1.0e-300)
        q13 = max((5.0 * q12 - val.bst * q11) / 4.0, 1.0e-300)
        # For charged pairs, Q14/Q15 can be carried from Debye-Huckel tables.
        # For electron-neutral pairs, keep the Magin-style fallback Q14=Q15=Q13.
        q14 = max(val.q14_m2, 1.0e-300) if val.q14_m2 is not None else q13
        q15 = max(val.q15_m2, 1.0e-300) if val.q15_m2 is not None else q13

        q11_eh.append(q11)
        q12_eh.append(q12)
        q13_eh.append(q13)
        q14_eh.append(q14)
        q15_eh.append(q15)

        w = x[sp] * q11
        electron_pair_weights_raw[f"e--{sp}"] = w
        wsum_e += w
        a_e_acc += w * val.ast
        b_e_acc += w * val.bst
        c_e_acc += w * val.cst
    electron_pair_weights = normalize_weights(electron_pair_weights_raw)

    a_e = a_e_acc / wsum_e if wsum_e > 0.0 else 1.0
    b_e = b_e_acc / wsum_e if wsum_e > 0.0 else 1.0
    c_e = c_e_acc / wsum_e if wsum_e > 0.0 else 1.0

    q22ee, q23ee, q24ee, _, _ = debye_huckel_repulsive_ee_integrals(t, n["e-"])
    x_e = max(x["e-"], 0.0)

    dot_l00 = sum(x[sp] * q11_eh[i] for i, sp in enumerate(HEAVY))
    dot_l01 = sum(x[sp] * (2.5 * q11_eh[i] - 3.0 * q12_eh[i]) for i, sp in enumerate(HEAVY))
    dot_l11 = sum(
        x[sp] * (6.25 * q11_eh[i] - 15.0 * q12_eh[i] + 12.0 * q13_eh[i]) for i, sp in enumerate(HEAVY)
    )
    dot_l02 = sum(
        x[sp] * (35.0 / 8.0 * q11_eh[i] - 10.5 * q12_eh[i] + 6.0 * q13_eh[i]) for i, sp in enumerate(HEAVY)
    )
    dot_l12 = sum(
        x[sp]
        * (175.0 / 16.0 * q11_eh[i] - 315.0 / 8.0 * q12_eh[i] + 57.0 * q13_eh[i] - 30.0 * q14_eh[i])
        for i, sp in enumerate(HEAVY)
    )
    dot_l22 = sum(
        x[sp]
        * (
            1225.0 / 64.0 * q11_eh[i]
            - 735.0 / 8.0 * q12_eh[i]
            + 199.5 * q13_eh[i]
            - 210.0 * q14_eh[i]
            + 90.0 * q15_eh[i]
        )
        for i, sp in enumerate(HEAVY)
    )

    lee3 = [[0.0 for _ in range(3)] for _ in range(3)]
    lee3[0][0] = dot_l00
    lee3[0][1] = dot_l01
    lee3[1][0] = dot_l01
    lee3[1][1] = dot_l11 + x_e * SQRT2 * q22ee
    lee3[0][2] = dot_l02
    lee3[2][0] = dot_l02
    lee3[1][2] = dot_l12 + x_e * SQRT2 * (1.75 * q22ee - 2.0 * q23ee)
    lee3[2][1] = lee3[1][2]
    lee3[2][2] = dot_l22 + x_e * SQRT2 * (77.0 / 16.0 * q22ee - 7.0 * q23ee + 5.0 * q24ee)

    lee_matrix_stat = matrix_metrics(lee3)
    p_pa = p * 101325.0
    leefac = 16.0 * p_pa / (3.0 * K_B * max(t, 1.0e-12)) * math.sqrt(MASS_E / (2.0 * PI * K_B * max(t, 1.0e-12)))
    de_sol = solve_linear_system(lee3, [1.0, 0.0, 0.0])
    dee3 = max(de_sol[0] / max(leefac, 1.0e-300), 0.0)
    sigma = max(n["e-"] * E_CHARGE * E_CHARGE / max(K_B * t, 1.0e-300) * dee3, 0.0)

    den_ke = lee3[1][1] * lee3[2][2] - lee3[1][2] * lee3[1][2]
    fac_ke = 75.0 * K_B / 64.0 * math.sqrt(2.0 * PI * K_B * t / MASS_E)
    if den_ke > 1.0e-300:
        k_e = fac_ke * x_e * lee3[2][2] / den_ke
    elif lee3[1][1] > 1.0e-300:
        k_e = fac_ke * x_e / lee3[1][1]
    else:
        k_e = 0.0
    k_e = max(k_e, 0.0)
    nu_e = n["e-"] * E_CHARGE * E_CHARGE / (MASS_E * sigma) if sigma > 0.0 else 0.0

    # --- Internal conductivity (Hirschfelder-Eucken) ---
    cp_frozen = 0.0
    h_species_mass: dict[str, float] = {}
    for sp in SPECIES:
        cp_i = interp1(t_pf, cp_mass_pf[sp], t)
        h_i = interp1(t_pf, h_mass_pf[sp], t)
        cp_frozen += y[sp] * cp_i
        h_species_mass[sp] = h_i

    cp_trans_total = (2.5 * K_B * n_total) / rho if n_total > 0.0 else 0.0
    cp_internal = max(cp_frozen - cp_trans_total, 0.0)
    k_int = mu_heavy * cp_internal / max(pr_int, 1.0e-8)

    # --- Reaction conductivity (Butler-Brokaw-like) ---
    diffusivity: dict[str, float] = {}
    for si in HEAVY:
        nu_i = 0.0
        for sj in SPECIES:
            val = db.get(si, sj, t, p)
            k11 = collision_rate_coeff(t, val.q11_m2, MASS[si], MASS[sj])
            nu_i += n[sj] * k11
        nu_i = max(nu_i, 1e-30)
        diffusivity[si] = K_B * t / (MASS[si] * nu_i)

    k_reac_raw = 0.0
    for sp in HEAVY:
        h_tilde = h_species_mass[sp] - h_mix
        k_reac_raw += diffusivity[sp] * h_tilde * y_deriv[sp]

    # Mass-fraction gradient form:
    #   k_reac = rho * sum_i(D_i * (h_i-h_mix) * dY_i/dT)
    k_reac_grad_mass = rho * k_reac_raw

    # Mole-fraction gradient form:
    #   k_reac = n_tot * sum_i(D_i * (h_i-h_mix)_particle * dx_i/dT)
    h_part = {sp: h_species_mass[sp] * MASS[sp] for sp in SPECIES}
    h_mix_part = sum(x[sp] * h_part[sp] for sp in SPECIES)
    k_reac_mole_sum = 0.0
    for sp in HEAVY:
        h_tilde_part = h_part[sp] - h_mix_part
        k_reac_mole_sum += diffusivity[sp] * h_tilde_part * x_deriv[sp]
    k_reac_grad_mole = n_total * k_reac_mole_sum

    cp_reac = cp_mix - cp_frozen
    k_reac_cp = mu_heavy * max(cp_reac, 0.0) / max(pr_heavy, 1e-8)

    if k_reac_model == "legacy_composite":
        k_reac_raw_model = k_reac_grad_mass + k_reac_cp
    elif k_reac_model == "butler_mass":
        k_reac_raw_model = k_reac_grad_mass
    elif k_reac_model == "butler_mole":
        k_reac_raw_model = k_reac_grad_mole
    else:  # butler_mole_cp
        k_reac_raw_model = k_reac_grad_mole + k_reac_cp

    k_reac = max(reaction_scale * k_reac_raw_model, 0.0)

    kappa = k_h + k_e + k_int + k_reac

    row = {
        "T_K": t,
        "P_atm": p,
        "mu_Pa_s": mu_heavy,
        "kappa_W_mK": kappa,
        "kappa_H_W_mK": k_h,
        "kappa_e_W_mK": k_e,
        "kappa_int_W_mK": k_int,
        "kappa_reac_W_mK": k_reac,
        "sigma_S_m": sigma,
        "nu_e_1_s": nu_e,
        "Astar_e": a_e,
        "Bstar_e": b_e,
        "Cstar_e": c_e,
        "cp_frozen_J_kgK": cp_frozen,
        "cp_reac_J_kgK": cp_reac,
        "kappa_reac_grad_mass_W_mK": k_reac_grad_mass,
        "kappa_reac_grad_mole_W_mK": k_reac_grad_mole,
        "kappa_reac_cp_W_mK": k_reac_cp,
    }

    heavy_pair_top, heavy_pair_top_share = dominant_pair(heavy_pair_weights)
    electron_pair_top, electron_pair_top_share = dominant_pair(electron_pair_weights)
    kappa_abs_sum = max(abs(k_h) + abs(k_e) + abs(k_int) + abs(k_reac), 1.0e-300)
    diagnostics: dict[str, object] = {
        "T_K": t,
        "P_atm": p,
        "collision_stage": {
            "heavy_pair_weight_share": heavy_pair_weights,
            "heavy_pair_top": heavy_pair_top,
            "heavy_pair_top_share": heavy_pair_top_share,
            "electron_pair_weight_share": electron_pair_weights,
            "electron_pair_top": electron_pair_top,
            "electron_pair_top_share": electron_pair_top_share,
        },
        "matrix_stage": {
            "mu": mu_matrix_stat,
            "k_H": kh_matrix_stat,
            "lee3": lee_matrix_stat,
        },
        "output_stage": {
            "mu_Pa_s": mu_heavy,
            "kappa_W_mK": kappa,
            "sigma_S_m": sigma,
            "kappa_share_H": k_h / kappa_abs_sum,
            "kappa_share_e": k_e / kappa_abs_sum,
            "kappa_share_int": k_int / kappa_abs_sum,
            "kappa_share_reac": k_reac / kappa_abs_sum,
        },
    }
    return row, diagnostics


def compute_composition_derivatives(
    rows: list[dict[str, float]],
) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    temps = [r["T_K"] for r in rows]

    y_by_sp: dict[str, list[float]] = {sp: [0.0] * len(rows) for sp in SPECIES}
    x_by_sp: dict[str, list[float]] = {sp: [0.0] * len(rows) for sp in SPECIES}
    for i, row in enumerate(rows):
        rho = 0.0
        n_total = 0.0
        for sp in SPECIES:
            rho += row[sp] * MASS[sp]
            n_total += row[sp]
        rho = max(rho, 1e-300)
        n_total = max(n_total, 1e-300)
        for sp in SPECIES:
            y_by_sp[sp][i] = row[sp] * MASS[sp] / rho
            x_by_sp[sp][i] = row[sp] / n_total

    dydt_by_sp = {sp: numeric_derivative(temps, y_by_sp[sp]) for sp in SPECIES}
    dxdt_by_sp = {sp: numeric_derivative(temps, x_by_sp[sp]) for sp in SPECIES}
    out_y: list[dict[str, float]] = []
    out_x: list[dict[str, float]] = []
    for i in range(len(rows)):
        out_y.append({sp: dydt_by_sp[sp][i] for sp in SPECIES})
        out_x.append({sp: dxdt_by_sp[sp][i] for sp in SPECIES})
    return out_y, out_x


def aggregate_state_diagnostics(diag_rows: list[dict[str, object]]) -> dict[str, object]:
    if not diag_rows:
        return {}

    def mean(values: list[float]) -> float:
        if not values:
            return 0.0
        return sum(values) / len(values)

    heavy_acc: dict[str, float] = defaultdict(float)
    elec_acc: dict[str, float] = defaultdict(float)
    mu_diag_ratio: list[float] = []
    mu_offdiag_to_diag: list[float] = []
    kh_diag_ratio: list[float] = []
    kh_offdiag_to_diag: list[float] = []
    lee_diag_ratio: list[float] = []
    lee_offdiag_to_diag: list[float] = []
    share_h: list[float] = []
    share_e: list[float] = []
    share_int: list[float] = []
    share_reac: list[float] = []

    for d in diag_rows:
        cstage = d["collision_stage"]
        for k, v in cstage["heavy_pair_weight_share"].items():
            heavy_acc[k] += float(v)
        for k, v in cstage["electron_pair_weight_share"].items():
            elec_acc[k] += float(v)

        mstage = d["matrix_stage"]
        mu_diag_ratio.append(float(mstage["mu"]["diag_ratio"]))
        mu_offdiag_to_diag.append(float(mstage["mu"]["offdiag_to_diag"]))
        kh_diag_ratio.append(float(mstage["k_H"]["diag_ratio"]))
        kh_offdiag_to_diag.append(float(mstage["k_H"]["offdiag_to_diag"]))
        lee_diag_ratio.append(float(mstage["lee3"]["diag_ratio"]))
        lee_offdiag_to_diag.append(float(mstage["lee3"]["offdiag_to_diag"]))

        ostage = d["output_stage"]
        share_h.append(float(ostage["kappa_share_H"]))
        share_e.append(float(ostage["kappa_share_e"]))
        share_int.append(float(ostage["kappa_share_int"]))
        share_reac.append(float(ostage["kappa_share_reac"]))

    n = float(len(diag_rows))
    heavy_mean = {k: v / n for k, v in heavy_acc.items()}
    elec_mean = {k: v / n for k, v in elec_acc.items()}
    heavy_sorted = sorted(heavy_mean.items(), key=lambda kv: kv[1], reverse=True)
    elec_sorted = sorted(elec_mean.items(), key=lambda kv: kv[1], reverse=True)

    return {
        "n_states": len(diag_rows),
        "collision_stage": {
            "heavy_pair_mean_share": heavy_mean,
            "electron_pair_mean_share": elec_mean,
            "heavy_pair_top3": heavy_sorted[:3],
            "electron_pair_top3": elec_sorted[:3],
        },
        "matrix_stage": {
            "mu_diag_ratio_mean": mean(mu_diag_ratio),
            "mu_offdiag_to_diag_mean": mean(mu_offdiag_to_diag),
            "kH_diag_ratio_mean": mean(kh_diag_ratio),
            "kH_offdiag_to_diag_mean": mean(kh_offdiag_to_diag),
            "lee3_diag_ratio_mean": mean(lee_diag_ratio),
            "lee3_offdiag_to_diag_mean": mean(lee_offdiag_to_diag),
        },
        "output_stage": {
            "kappa_share_H_mean": mean(share_h),
            "kappa_share_e_mean": mean(share_e),
            "kappa_share_int_mean": mean(share_int),
            "kappa_share_reac_mean": mean(share_reac),
        },
    }


def write_csv(path: Path, rows: list[dict[str, float]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in fieldnames})


def pressure_tag(p_atm: float) -> str:
    return f"{p_atm:g}".replace("-", "m").replace(".", "p")


def generate_plots(
    combined_csv: Path, plots_dir: Path, pressures: list[float], matf_csv: Path | None
) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    pressure_terms = []
    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    for idx, p in enumerate(pressures):
        if abs(p - 0.1) < 1e-12:
            label = "0.1 atm"
        elif abs(p - 1.0) < 1e-12:
            label = "1 atm"
        elif abs(p - 4.0) < 1e-12:
            label = "4 atm"
        else:
            label = f"{p:g} atm"
        pressure_terms.append((p, colors[idx % len(colors)], label))

    summary_scripts = [
        ("mu_vs_T.gnuplot", 3, 9, "Viscosity", "mu [Pa s]", plots_dir / "mu_vs_T.png"),
        (
            "kappa_vs_T.gnuplot",
            4,
            10,
            "Thermal Conductivity",
            "kappa [W/(m K)]",
            plots_dir / "kappa_vs_T.png",
        ),
        (
            "sigma_vs_T.gnuplot",
            9,
            11,
            "Electrical Conductivity",
            "sigma [S/m]",
            plots_dir / "sigma_vs_T.png",
        ),
    ]

    for script_name, col, matf_col, title, ylabel, png_path in summary_scripts:
        lines = []
        for p, color, label in pressure_terms:
            lines.append(
                f"'{combined_csv}' u 1:(abs($2-{p:.16g})<1e-12 ? ${col} : 1/0) "
                f"w l lw 2 lc rgb '{color}' title 'This work {label}'"
            )
            if matf_csv is not None and matf_csv.exists():
                lines.append(
                    f"'{matf_csv}' u 1:(abs($2-{p:.16g})<1e-12 ? ${matf_col} : 1/0) "
                    f"w l lw 2 dt 2 lc rgb '{color}' title 'MATF {label}'"
                )
        script = (
            "set datafile separator ','\n"
            "set terminal pngcairo size 1400,900 enhanced\n"
            f"set output '{png_path}'\n"
            "set grid\n"
            "set key outside right top\n"
            "set xlabel 'Temperature [K]'\n"
            f"set ylabel '{ylabel}'\n"
            f"set title '{title}'\n"
            "set format y '%.3e'\n"
            "plot \\\n"
            + ", \\\n".join(lines)
            + "\n"
        )
        path = plots_dir / script_name
        path.write_text(script, encoding="utf-8")
        subprocess.run(["gnuplot", str(path)], check=True)

    for p, _, label in pressure_terms:
        ptag = pressure_tag(p)
        png_path = plots_dir / f"kappa_components_{ptag}atm.png"
        script_path = plots_dir / f"kappa_components_{ptag}atm.gnuplot"
        script = (
            "set datafile separator ','\n"
            "set terminal pngcairo size 1500,950 enhanced\n"
            f"set output '{png_path}'\n"
            "set grid\n"
            "set key outside right top\n"
            "set xlabel 'Temperature [K]'\n"
            "set ylabel 'kappa component [W/(m K)]'\n"
            f"set title 'Thermal Conductivity Components ({label})'\n"
            "set format y '%.3e'\n"
            "plot \\\n"
            f"'{combined_csv}' u 1:(abs($2-{p:.16g})<1e-12 ? $4 : 1/0) w l lw 2 lc rgb '#000000' title 'kappa total', \\\n"
            f"'{combined_csv}' u 1:(abs($2-{p:.16g})<1e-12 ? $5 : 1/0) w l lw 2 lc rgb '#1f77b4' title 'k_H', \\\n"
            f"'{combined_csv}' u 1:(abs($2-{p:.16g})<1e-12 ? $6 : 1/0) w l lw 2 lc rgb '#2ca02c' title 'k_e', \\\n"
            f"'{combined_csv}' u 1:(abs($2-{p:.16g})<1e-12 ? $7 : 1/0) w l lw 2 lc rgb '#ff7f0e' title 'k_int', \\\n"
            f"'{combined_csv}' u 1:(abs($2-{p:.16g})<1e-12 ? $8 : 1/0) w l lw 2 lc rgb '#d62728' title 'k_reac'\n"
        )
        script_path.write_text(script, encoding="utf-8")
        subprocess.run(["gnuplot", str(script_path)], check=True)


def main() -> None:
    args = parse_args()
    apply_murphy_profile(args)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / "plots"

    equilibrium_by_p, _, pressures = load_equilibrium(args.equilibrium_csv.resolve())
    thermo_by_pt = load_thermo(args.thermo_csv.resolve())
    t_pf, h_mass_pf, cp_mass_pf, _ = load_species_frozen_thermo(args.partition_json.resolve())

    records, mapping = load_collision_records(args.collision_csv.resolve())
    if not records:
        raise RuntimeError("No usable collision records were loaded")
    db = CollisionDatabase(records)

    combined_rows: list[dict[str, float]] = []
    per_pressure: dict[float, list[dict[str, float]]] = {}
    diag_by_pressure: dict[float, list[dict[str, object]]] = defaultdict(list)

    for p in pressures:
        states = equilibrium_by_p[p]
        dydt_rows, dxdt_rows = compute_composition_derivatives(states)
        rows_p: list[dict[str, float]] = []

        for idx, state in enumerate(states):
            t = state["T_K"]
            thermo = thermo_by_pt.get((p, t))
            if thermo is None:
                raise KeyError(f"Missing thermo row for P={p}, T={t}")

            row, diag = compute_transport_for_state(
                state=state,
                thermo=thermo,
                y_deriv=dydt_rows[idx],
                x_deriv=dxdt_rows[idx],
                t_pf=t_pf,
                h_mass_pf=h_mass_pf,
                cp_mass_pf=cp_mass_pf,
                db=db,
                pr_heavy=args.pr_heavy,
                pr_int=args.pr_int,
                reaction_scale=args.reaction_scale,
                k_reac_model=args.k_reac_model,
            )
            rows_p.append(row)
            combined_rows.append(row)
            diag_by_pressure[p].append(diag)

        per_pressure[p] = rows_p

    combined_rows.sort(key=lambda r: (r["P_atm"], r["T_K"]))

    fields = [
        "T_K",
        "P_atm",
        "mu_Pa_s",
        "kappa_W_mK",
        "kappa_H_W_mK",
        "kappa_e_W_mK",
        "kappa_int_W_mK",
        "kappa_reac_W_mK",
        "sigma_S_m",
        "nu_e_1_s",
        "Astar_e",
        "Bstar_e",
        "Cstar_e",
        "cp_frozen_J_kgK",
        "cp_reac_J_kgK",
        "kappa_reac_grad_mass_W_mK",
        "kappa_reac_grad_mole_W_mK",
        "kappa_reac_cp_W_mK",
    ]

    combined_csv = output_dir / "argon_transport_0p1_1_4atm.csv"
    write_csv(combined_csv, combined_rows, fields)

    for p in pressures:
        ptag = pressure_tag(p)
        per_csv = output_dir / f"argon_transport_{ptag}atm.csv"
        write_csv(per_csv, per_pressure[p], fields)

    matf_reference_csv = args.matf_reference_csv.resolve()
    if not matf_reference_csv.exists():
        print(f"[WARN] MATF reference CSV not found: {matf_reference_csv}")
        matf_reference_csv = None

    sensitivity_report: dict[str, object] = {}
    for p in pressures:
        key = f"{p:g} atm"
        sensitivity_report[key] = aggregate_state_diagnostics(diag_by_pressure[p])

    sensitivity_path = output_dir / "argon_transport_sensitivity_report.json"
    sensitivity_path.write_text(
        json.dumps(sensitivity_report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    for p in pressures:
        key = f"{p:g} atm"
        srep = sensitivity_report.get(key, {})
        cstage = srep.get("collision_stage", {})
        ostage = srep.get("output_stage", {})
        top_heavy = cstage.get("heavy_pair_top3", [])
        if top_heavy:
            pair0, share0 = top_heavy[0]
        else:
            pair0, share0 = "", 0.0
        print(
            "[SENS] "
            f"{key}: heavy-top={pair0} ({100.0*float(share0):.2f}%), "
            f"kappa shares [H={100.0*float(ostage.get('kappa_share_H_mean', 0.0)):.2f}%, "
            f"e={100.0*float(ostage.get('kappa_share_e_mean', 0.0)):.2f}%, "
            f"int={100.0*float(ostage.get('kappa_share_int_mean', 0.0)):.2f}%, "
            f"reac={100.0*float(ostage.get('kappa_share_reac_mean', 0.0)):.2f}%]"
        )

    peak_summary: dict[str, dict[str, float]] = {}
    for p in pressures:
        rows = per_pressure[p]
        peak_mu = max(rows, key=lambda r: r["mu_Pa_s"])
        peak_kr = max(rows, key=lambda r: r["kappa_reac_W_mK"])
        peak_k = max(rows, key=lambda r: r["kappa_W_mK"])
        peak_sigma = max(rows, key=lambda r: r["sigma_S_m"])
        peak_summary[f"{p:g} atm"] = {
            "mu_peak_T_K": peak_mu["T_K"],
            "mu_peak_Pa_s": peak_mu["mu_Pa_s"],
            "kappa_peak_T_K": peak_k["T_K"],
            "kappa_peak_W_mK": peak_k["kappa_W_mK"],
            "kappa_reac_peak_T_K": peak_kr["T_K"],
            "kappa_reac_peak_W_mK": peak_kr["kappa_reac_W_mK"],
            "sigma_peak_T_K": peak_sigma["T_K"],
            "sigma_peak_S_m": peak_sigma["sigma_S_m"],
        }

    metadata = {
        "inputs": {
            "equilibrium_csv": str(args.equilibrium_csv.resolve()),
            "thermo_csv": str(args.thermo_csv.resolve()),
            "collision_csv": str(args.collision_csv.resolve()),
            "partition_json": str(args.partition_json.resolve()),
            "matf_reference_csv": (
                str(args.matf_reference_csv.resolve())
                if args.matf_reference_csv is not None
                else ""
            ),
        },
        "collision_mapping": mapping,
        "collision_rows_loaded": len(records),
        "fallback_counts": dict(sorted(db.fallback_counts.items())),
        "approximations": {
            "viscosity": "Devoto 1st-order Chapman-Enskog heavy-species coefficient matrix",
            "kappa_H": "Devoto/Muckenfuss-Curtiss 2nd-order heavy-species coefficient matrix",
            "kappa_e": "Devoto 3rd-order electron Lee matrix",
            "kappa_int": "Hirschfelder-Eucken internal conductivity (algebraic)",
            "kappa_reac": "Configurable Butler-Brokaw-style closure selected by --k-reac-model",
            "sigma": "Devoto 3rd-order electron conductivity from Lee matrix inverse",
            "electron_missing_integrals": "Q14ei/Q15ei read from Phase-4 CSV when available (including e-Ar regularized closure rows); fallback Q14ei=Q15ei=Q13ei is used only when absent.",
        },
        "constants": {
            "reaction_scale": args.reaction_scale,
            "pr_heavy": args.pr_heavy,
            "pr_int": args.pr_int,
            "k_reac_model": args.k_reac_model,
            "murphy_strict": args.murphy_strict,
        },
        "peaks": peak_summary,
        "notes": [
            "kappa_reac peak should be inspected in ionization region (typically O(1e4-2e4 K), pressure dependent).",
            "compute_transport_for_state() uses strict Devoto-style coefficient matrices for mu, kappa_H, kappa_e, and sigma.",
        ],
        "outputs": {
            "combined_csv": str(combined_csv),
            "per_pressure_csv": [
                str(output_dir / f"argon_transport_{pressure_tag(p)}atm.csv") for p in pressures
            ],
            "sensitivity_report_json": str(sensitivity_path),
            "plots_dir": str(plots_dir),
            "matf_overlay_used": matf_reference_csv is not None,
        },
    }
    meta_path = output_dir / "argon_transport_metadata.json"
    meta_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    if not args.skip_plots:
        try:
            generate_plots(
                combined_csv=combined_csv,
                plots_dir=plots_dir,
                pressures=pressures,
                matf_csv=matf_reference_csv,
            )
        except FileNotFoundError:
            print("[WARN] gnuplot not found. Skipping plots.")
        except subprocess.CalledProcessError as exc:
            print(f"[WARN] gnuplot failed ({exc}). Skipping plots.")

    print(f"[OK] Wrote: {combined_csv}")
    for p in pressures:
        print(f"[OK] Wrote: {output_dir / f'argon_transport_{pressure_tag(p)}atm.csv'}")
    print(f"[OK] Wrote: {sensitivity_path}")
    print(f"[OK] Wrote: {meta_path}")
    if not args.skip_plots:
        print(f"[OK] Wrote plots in: {plots_dir}")


if __name__ == "__main__":
    main()
