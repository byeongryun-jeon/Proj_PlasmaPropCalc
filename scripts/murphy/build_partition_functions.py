#!/usr/bin/env python3
"""Build internal partition-function lookup tables for Ar plasma species.

This script expects NIST-level CSV files named Ar0.csv ... Ar4.csv, extracts
`g` and `Level (eV)` columns, and computes:

    Q(T) = sum_i g_i * exp(-epsilon_i / (k_B * T))

for each species. No ionization-potential-lowering cutoff is applied.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

try:
    import numpy as np  # Optional; used only for NPZ export.
except ModuleNotFoundError:
    np = None


KB_EV_PER_K = 8.617333262145e-5  # Boltzmann constant [eV/K]
NUM_PATTERN = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
CHARGE_TO_SPECIES = {
    0: "Ar",
    1: "Ar+",
    2: "Ar2+",
    3: "Ar3+",
    4: "Ar4+",
}


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Compute Ar/Ar+/Ar2+/Ar3+/Ar4+/e- partition-function tables."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_root / "data" / "raw" / "murphy" / "nist_levels" / "argon",
        help="Directory containing Ar0.csv ... Ar4.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_root / "data" / "processed" / "partition_functions",
        help="Directory to write lookup-table outputs",
    )
    parser.add_argument("--t-min", type=float, default=300.0, help="Minimum T [K]")
    parser.add_argument("--t-max", type=float, default=30000.0, help="Maximum T [K]")
    parser.add_argument("--t-step", type=float, default=200.0, help="Temperature step [K]")
    parser.add_argument(
        "--exclude-t-max",
        action="store_true",
        help="Do not force-include t-max if it is off-grid by t-step.",
    )
    return parser.parse_args()


def extract_numeric(value: object) -> float:
    """Extract the first numeric token from a raw CSV cell string."""
    if value is None:
        raise ValueError("Empty value")
    text = str(value).strip()
    match = NUM_PATTERN.search(text)
    if not match:
        raise ValueError(f"No numeric token found in cell: {text!r}")
    return float(match.group(0))


def load_levels(
    csv_path: Path,
) -> tuple[list[float], list[float], dict[str, int], list[float]]:
    """Load g and level(eV) arrays from one species CSV file.

    Rows are skipped when:
    - g is missing
    - g is not numeric
    - level is missing
    - level is not numeric

    Additionally, rows with missing g but valid level are captured as
    ionization-limit candidates for this charge state.
    """
    g_values: list[float] = []
    eps_values: list[float] = []
    ionization_limit_candidates_eV: list[float] = []
    skipped = {
        "missing_g": 0,
        "invalid_g": 0,
        "missing_level_eV": 0,
        "invalid_level_eV": 0,
    }

    with csv_path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        if reader.fieldnames is None:
            raise ValueError(f"{csv_path} has no header row")
        for row in reader:
            raw_g = (row.get("g") or "").strip()
            raw_eps = (row.get("Level (eV)") or "").strip()

            if not raw_g:
                if raw_eps:
                    try:
                        ionization_limit_candidates_eV.append(extract_numeric(raw_eps))
                    except ValueError:
                        pass
                skipped["missing_g"] += 1
                continue
            if not raw_eps:
                skipped["missing_level_eV"] += 1
                continue

            try:
                g = extract_numeric(raw_g)
            except ValueError:
                skipped["invalid_g"] += 1
                continue

            try:
                eps = extract_numeric(raw_eps)
            except ValueError:
                skipped["invalid_level_eV"] += 1
                continue
            g_values.append(g)
            eps_values.append(eps)

    if not g_values:
        raise ValueError(f"{csv_path} contains no usable data rows")
    return g_values, eps_values, skipped, ionization_limit_candidates_eV


def compute_q_of_t(g: list[float], eps_ev: list[float], temperatures: list[float]) -> list[float]:
    """Compute Q(T) for one species."""
    q_values: list[float] = []
    for temp in temperatures:
        denom = KB_EV_PER_K * temp
        total = 0.0
        for gi, ei in zip(g, eps_ev):
            total += gi * math.exp(-ei / denom)
        q_values.append(total)
    return q_values


def build_temperature_grid(
    t_min: float, t_max: float, t_step: float, include_t_max: bool
) -> list[float]:
    if t_step <= 0.0:
        raise ValueError("t-step must be > 0")
    if t_min <= 0.0:
        raise ValueError("t-min must be > 0 K")
    if t_max < t_min:
        raise ValueError("t-max must be >= t-min")

    temps: list[float] = []
    t = t_min
    while t <= t_max + 1e-12:
        temps.append(float(t))
        t += t_step

    if include_t_max and (not temps or abs(temps[-1] - t_max) > 1e-9):
        temps.append(float(t_max))

    return temps


def charge_from_filename(csv_path: Path) -> int:
    match = re.fullmatch(r"Ar([0-4])\.csv", csv_path.name)
    if not match:
        raise ValueError(f"Unexpected file name: {csv_path.name}")
    return int(match.group(1))


def key_from_species(species: str) -> str:
    if species == "e-":
        return "e"
    return species.replace("+", "_p")


def write_csv(
    csv_path: Path, temperatures: list[float], species_names: list[str], q_by_species: dict[str, list[float]]
) -> None:
    header = ["T_K"] + [f"Q_{key_from_species(name)}" for name in species_names]
    with csv_path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(header)
        for idx, temp in enumerate(temperatures):
            row = [f"{temp:.6f}"]
            for species in species_names:
                row.append(f"{q_by_species[species][idx]:.16e}")
            writer.writerow(row)


def maybe_write_npz(
    npz_path: Path, temperatures: list[float], species_names: list[str], q_by_species: dict[str, list[float]]
) -> bool:
    if np is None:
        return False

    q_matrix = np.vstack([np.asarray(q_by_species[name], dtype=np.float64) for name in species_names])
    np.savez_compressed(
        npz_path,
        temperature_K=np.asarray(temperatures, dtype=np.float64),
        species_names=np.asarray(species_names),
        Q_matrix=q_matrix,
        kb_eV_per_K=np.float64(KB_EV_PER_K),
    )
    return True


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_files = sorted(input_dir.glob("Ar[0-4].csv"), key=charge_from_filename)
    if len(csv_files) != 5:
        raise FileNotFoundError(
            f"Expected 5 files Ar0.csv..Ar4.csv in {input_dir}, found {len(csv_files)}"
        )

    temperatures = build_temperature_grid(
        t_min=args.t_min,
        t_max=args.t_max,
        t_step=args.t_step,
        include_t_max=not args.exclude_t_max,
    )

    species_names: list[str] = []
    q_by_species: dict[str, list[float]] = {}
    source_summary: dict[str, dict[str, object]] = {}

    stage_ionization_eV_by_charge: dict[int, float] = {}

    for csv_file in csv_files:
        charge = charge_from_filename(csv_file)
        species = CHARGE_TO_SPECIES[charge]
        g, eps, skipped, ionization_limit_candidates_eV = load_levels(csv_file)
        q_t = compute_q_of_t(g=g, eps_ev=eps, temperatures=temperatures)

        species_names.append(species)
        q_by_species[species] = q_t
        skipped_rows = int(sum(skipped.values()))

        selected_stage_ionization_eV = None
        if ionization_limit_candidates_eV:
            selected_stage_ionization_eV = ionization_limit_candidates_eV[0]
            stage_ionization_eV_by_charge[charge] = selected_stage_ionization_eV

        source_summary[species] = {
            "source_file": str(csv_file),
            "n_levels": len(eps),
            "skipped_rows": skipped_rows,
            "skipped_breakdown": skipped,
            "ionization_limit_candidates_eV": ionization_limit_candidates_eV,
            "selected_stage_ionization_eV": selected_stage_ionization_eV,
            "min_level_eV": min(eps),
            "max_level_eV": max(eps),
        }

    cumulative_ionization_eV_by_species: dict[str, float] = {"Ar": 0.0}
    running_eion = 0.0
    for charge in (0, 1, 2, 3):
        if charge not in stage_ionization_eV_by_charge:
            raise ValueError(
                "Missing stage ionization limit in raw data for "
                f"{CHARGE_TO_SPECIES[charge]} (expected g-missing row with level)."
            )
        running_eion += stage_ionization_eV_by_charge[charge]
        cumulative_ionization_eV_by_species[CHARGE_TO_SPECIES[charge + 1]] = running_eion
    cumulative_ionization_eV_by_species["e-"] = 0.0

    # Electron internal partition function (temperature-independent)
    species_names.append("e-")
    q_by_species["e-"] = [2.0] * len(temperatures)
    source_summary["e-"] = {
        "source_file": "constant",
        "n_levels": 1,
        "min_level_eV": 0.0,
        "max_level_eV": 0.0,
        "Q_constant": 2.0,
    }

    csv_path = output_dir / "partition_functions_Ar_system.csv"
    json_path = output_dir / "partition_functions_Ar_system.json"
    npz_path = output_dir / "partition_functions_Ar_system.npz"
    meta_path = output_dir / "partition_functions_Ar_system_metadata.json"

    write_csv(csv_path, temperatures, species_names, q_by_species)
    wrote_npz = maybe_write_npz(npz_path, temperatures, species_names, q_by_species)

    table_payload = {
        "temperature_K": temperatures,
        "species_order": species_names,
        "Q": {species: q_by_species[species] for species in species_names},
        "stage_ionization_eV_by_charge": stage_ionization_eV_by_charge,
        "cumulative_ionization_eV_by_species": cumulative_ionization_eV_by_species,
        "kb_eV_per_K": KB_EV_PER_K,
        "formula": "Q(T) = sum_i g_i * exp(-epsilon_i / (k_B * T))",
    }
    json_path.write_text(json.dumps(table_payload, ensure_ascii=False), encoding="utf-8")

    metadata = {
        "formula": "Q(T) = sum_i g_i * exp(-epsilon_i / (k_B * T))",
        "kb_eV_per_K": KB_EV_PER_K,
        "temperature_grid_K": {
            "t_min": temperatures[0],
            "t_max": temperatures[-1],
            "n_points": len(temperatures),
            "nominal_step": args.t_step,
            "forced_include_t_max": not args.exclude_t_max,
        },
        "species_order": species_names,
        "stage_ionization_eV_by_charge": stage_ionization_eV_by_charge,
        "cumulative_ionization_eV_by_species": cumulative_ionization_eV_by_species,
        "notes": [
            "All energy levels in each CSV are included (no ionization potential lowering cutoff).",
            "Designed for operating pressures <= 4 atm per project instruction.",
        ],
        "sources": source_summary,
        "outputs": {
            "csv": str(csv_path),
            "json": str(json_path),
            "npz": str(npz_path) if wrote_npz else None,
        },
        "numpy_available": wrote_npz,
    }
    meta_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"[OK] Wrote: {csv_path}")
    print(f"[OK] Wrote: {json_path}")
    print(f"[OK] Wrote: {meta_path}")
    if wrote_npz:
        print(f"[OK] Wrote: {npz_path}")
    else:
        print("[INFO] NumPy not available, skipped NPZ output.")
    print(f"[INFO] Species: {', '.join(species_names)}")
    print(
        "[INFO] Temperature range: "
        f"{temperatures[0]:.1f} K .. {temperatures[-1]:.1f} K ({len(temperatures)} points)"
    )


if __name__ == "__main__":
    main()
