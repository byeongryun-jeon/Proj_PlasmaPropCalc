#!/usr/bin/env python3
"""Build thermodynamic lookup tables (rho, h, Cp) for Ar LTE plasma.

Inputs:
    - partition_functions_Ar_system.json  (Q_int, cumulative ionization energies)
    - argon_lte_equilibrium_0p1_1_4atm.csv (equilibrium number densities)

Outputs:
    - per-pressure CSV tables for rho [kg/m^3], h [J/kg], Cp [J/(kg K)]
    - combined CSV table
    - metadata JSON summary
    - verification plots (rho/h/Cp vs T) as PNG via gnuplot
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import mimetypes
import os
import shutil
import subprocess
from pathlib import Path


# --- Constants ---
KB_EV_PER_K = 8.617333262e-5
EV_TO_J = 1.602176634e-19
MASS_AR = 6.6335209e-26
MASS_E = 9.1093837e-31

SPECIES = ["Ar", "Ar+", "Ar2+", "Ar3+", "Ar4+", "e-"]
EQ_COL_TO_SPECIES = {
    "n_Ar": "Ar",
    "n_Ar_p": "Ar+",
    "n_Ar2_p": "Ar2+",
    "n_Ar3_p": "Ar3+",
    "n_Ar4_p": "Ar4+",
    "n_e": "e-",
}
CP_PEAK_TMIN = 0.0
CP_PEAK_TMAX = 30000.0


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Compute rho/h/Cp lookup tables from Phase 1/2 outputs."
    )
    parser.add_argument(
        "--partition-json",
        type=Path,
        default=root / "data" / "processed" / "partition_functions" / "partition_functions_Ar_system.json",
        help="Partition-function JSON path.",
    )
    parser.add_argument(
        "--equilibrium-csv",
        type=Path,
        default=root / "data" / "processed" / "equilibrium" / "argon_lte_equilibrium_0p1_1_4atm.csv",
        help="Equilibrium composition CSV path.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "data" / "processed" / "thermo",
        help="Output directory for thermo tables and plots.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip gnuplot figure generation.",
    )
    parser.add_argument(
        "--upload-to-gdrive",
        action="store_true",
        help="Upload generated CSV/PNG outputs to Google Drive.",
    )
    parser.add_argument(
        "--gdrive-credentials",
        type=Path,
        default=root / "credentials.json",
        help="Path to Google credentials JSON (service account or OAuth client).",
    )
    parser.add_argument(
        "--gdrive-folder-id",
        type=str,
        default=None,
        help="Google Drive folder ID. If omitted, uses GDRIVE_FOLDER_ID env var.",
    )
    parser.add_argument(
        "--gdrive-token",
        type=Path,
        default=root / ".gdrive_token.json",
        help="Token cache path for OAuth credentials (ignored for service accounts).",
    )
    parser.add_argument(
        "--gdrive-subfolder-name",
        type=str,
        default=None,
        help="Optional subfolder name to create/use under the target Google Drive folder.",
    )
    return parser.parse_args()


def numerical_derivative(x: list[float], y: list[float]) -> list[float]:
    """Finite-difference derivative dy/dx with central difference interior."""
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")
    n = len(x)
    if n < 2:
        raise ValueError("Need at least 2 points for derivative.")

    dy = [0.0] * n
    dx0 = x[1] - x[0]
    dxn = x[-1] - x[-2]
    if dx0 == 0.0 or dxn == 0.0:
        raise ValueError("Repeated x values are not allowed.")
    dy[0] = (y[1] - y[0]) / dx0
    dy[-1] = (y[-1] - y[-2]) / dxn

    for i in range(1, n - 1):
        dx = x[i + 1] - x[i - 1]
        if dx == 0.0:
            raise ValueError("Repeated x values are not allowed.")
        dy[i] = (y[i + 1] - y[i - 1]) / dx
    return dy


def interp_linear(x_grid: list[float], y_grid: list[float], x: float) -> float:
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


def pressure_tag(p_atm: float) -> str:
    text = f"{p_atm:g}"
    return text.replace("-", "m").replace(".", "p")


def load_partition_data(partition_json: Path) -> tuple[list[float], dict[str, list[float]], dict[str, float]]:
    payload = json.loads(partition_json.read_text(encoding="utf-8"))
    temperature_grid = [float(t) for t in payload["temperature_K"]]
    q = payload["Q"]
    q_by_species = {sp: [float(v) for v in q[sp]] for sp in SPECIES}
    eion = payload["cumulative_ionization_eV_by_species"]
    eion_by_species = {sp: float(eion[sp]) for sp in SPECIES}
    return temperature_grid, q_by_species, eion_by_species


def load_equilibrium_data(
    equilibrium_csv: Path,
) -> tuple[list[float], dict[float, list[dict[str, float]]], list[float]]:
    by_pressure: dict[float, list[dict[str, float]]] = {}
    pressures: set[float] = set()
    temp_set: set[float] = set()

    with equilibrium_csv.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            t = float(row["T_K"])
            p = float(row["P_atm"])
            temp_set.add(t)
            pressures.add(p)

            n_species = {
                EQ_COL_TO_SPECIES[col]: float(row[col]) for col in EQ_COL_TO_SPECIES
            }
            by_pressure.setdefault(p, []).append({"T_K": t, **n_species})

    if not by_pressure:
        raise ValueError(f"No rows found in {equilibrium_csv}")

    for p in by_pressure:
        by_pressure[p].sort(key=lambda r: r["T_K"])
    temperatures = sorted(temp_set)
    pressure_order = sorted(pressures)
    return temperatures, by_pressure, pressure_order


def compute_species_enthalpy_j_per_particle(
    temperature_pf: list[float],
    q_by_species: dict[str, list[float]],
    eion_by_species_eV: dict[str, float],
) -> tuple[dict[str, list[float]], dict[str, list[float]]]:
    """Return H_i(T) [J/particle] and dlnQ/dT arrays on partition temperature grid."""
    dlnq_dt: dict[str, list[float]] = {}
    h_species_j: dict[str, list[float]] = {}

    for sp in SPECIES:
        q_vals = q_by_species[sp]
        lnq = [math.log(max(v, 1e-300)) for v in q_vals]
        dlnq = numerical_derivative(temperature_pf, lnq)
        dlnq_dt[sp] = dlnq

        h_ev = [
            2.5 * KB_EV_PER_K * t + KB_EV_PER_K * (t**2) * dlnq_i + eion_by_species_eV[sp]
            for t, dlnq_i in zip(temperature_pf, dlnq)
        ]
        h_species_j[sp] = [val * EV_TO_J for val in h_ev]

    return h_species_j, dlnq_dt


def build_h_lookup_for_temperature(
    temperature_pf: list[float], h_species_j: dict[str, list[float]], temperature: float
) -> dict[str, float]:
    return {
        sp: interp_linear(temperature_pf, h_species_j[sp], temperature) for sp in SPECIES
    }


def cp_peaks_between(
    rows: list[dict[str, float]], t_min: float, t_max: float
) -> tuple[float, float]:
    sub = [r for r in rows if t_min <= r["T_K"] <= t_max]
    if not sub:
        return float("nan"), float("nan")
    peak = max(sub, key=lambda r: r["cp_J_kgK"])
    return peak["T_K"], peak["cp_J_kgK"]


def write_csv(path: Path, rows: list[dict[str, float]], include_pressure: bool) -> None:
    if include_pressure:
        fields = ["T_K", "P_atm", "rho_kg_m3", "h_J_kg", "cp_J_kgK"]
    else:
        fields = ["T_K", "rho_kg_m3", "h_J_kg", "cp_J_kgK"]
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in fields})


def write_cp_peak_csv(path: Path, peak_rows: list[dict[str, float]]) -> None:
    fields = ["P_atm", "cp_peak_T_K_0_30000", "cp_peak_J_kgK_0_30000"]
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        writer.writerows(peak_rows)


def _drive_imports() -> tuple[object, object, object, object, object, object, object]:
    """Lazy import Google Drive modules so non-upload runs keep minimal deps."""
    try:
        from googleapiclient.discovery import build  # type: ignore
        from googleapiclient.errors import HttpError  # type: ignore
        from googleapiclient.http import MediaFileUpload  # type: ignore
        from google.auth.transport.requests import Request  # type: ignore
        from google.oauth2 import service_account  # type: ignore
        from google.oauth2.credentials import Credentials as UserCredentials  # type: ignore
        from google_auth_oauthlib.flow import InstalledAppFlow  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Google API libraries are not installed in this Python environment. "
            "Use the project venv: "
            "/home/brjeon/Proj_PlasmaPropCalc/.venv/bin/python"
        ) from exc
    return (
        build,
        HttpError,
        MediaFileUpload,
        Request,
        service_account,
        UserCredentials,
        InstalledAppFlow,
    )


def build_drive_service(credentials_path: Path, token_path: Path) -> object:
    build, _, _, Request, service_account, UserCredentials, InstalledAppFlow = _drive_imports()

    if not credentials_path.exists():
        raise FileNotFoundError(f"Google credentials file not found: {credentials_path}")

    scopes = ["https://www.googleapis.com/auth/drive"]
    payload = json.loads(credentials_path.read_text(encoding="utf-8"))

    creds = None
    if payload.get("type") == "service_account":
        creds = service_account.Credentials.from_service_account_file(
            str(credentials_path), scopes=scopes
        )
    elif "installed" in payload or "web" in payload:
        if token_path.exists():
            creds = UserCredentials.from_authorized_user_file(str(token_path), scopes=scopes)
        if creds is None or not creds.valid:
            if creds is not None and creds.expired and creds.refresh_token:
                creds.refresh(Request())
            else:
                flow = InstalledAppFlow.from_client_secrets_file(
                    str(credentials_path), scopes=scopes
                )
                creds = flow.run_local_server(port=0)
            token_path.write_text(creds.to_json(), encoding="utf-8")
    else:
        raise ValueError(
            "Unsupported Google credentials format. "
            "Expected service account JSON or OAuth client JSON."
        )

    return build("drive", "v3", credentials=creds, cache_discovery=False)


def _drive_query_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("'", "\\'")


def upload_file_to_drive(service: object, local_path: Path, folder_id: str) -> dict[str, str]:
    _, HttpError, MediaFileUpload, _, _, _, _ = _drive_imports()
    if not local_path.exists():
        raise FileNotFoundError(f"Upload target not found: {local_path}")

    mime_type, _ = mimetypes.guess_type(str(local_path))
    if mime_type is None:
        mime_type = "application/octet-stream"

    escaped_name = _drive_query_escape(local_path.name)
    query = f"'{folder_id}' in parents and name = '{escaped_name}' and trashed = false"
    listed = (
        service.files()
        .list(
            q=query,
            fields="files(id,name)",
            pageSize=1,
            supportsAllDrives=True,
            includeItemsFromAllDrives=True,
        )
        .execute()
    )
    existing = listed.get("files", [])
    media = MediaFileUpload(str(local_path), mimetype=mime_type, resumable=False)

    try:
        if existing:
            file_id = existing[0]["id"]
            resp = (
                service.files()
                .update(
                    fileId=file_id,
                    media_body=media,
                    fields="id,name,mimeType,webViewLink,webContentLink",
                    supportsAllDrives=True,
                )
                .execute()
            )
            action = "updated"
        else:
            resp = (
                service.files()
                .create(
                    body={"name": local_path.name, "parents": [folder_id]},
                    media_body=media,
                    fields="id,name,mimeType,webViewLink,webContentLink",
                    supportsAllDrives=True,
                )
                .execute()
            )
            action = "created"
    except HttpError as exc:
        err = str(exc)
        if "storageQuotaExceeded" in err:
            raise RuntimeError(
                "Google Drive 업로드 실패: 현재 credentials가 서비스 계정이며 "
                "대상 폴더가 Shared Drive가 아니어서 quota 제한이 발생했습니다. "
                "해결: (1) OAuth client credentials 사용 또는 "
                "(2) Shared Drive 폴더 ID 사용."
            ) from exc
        raise

    file_id = resp["id"]
    return {
        "action": action,
        "name": resp.get("name", local_path.name),
        "id": file_id,
        "mimeType": resp.get("mimeType", mime_type),
        "webViewLink": resp.get("webViewLink", f"https://drive.google.com/file/d/{file_id}/view"),
        "webContentLink": resp.get("webContentLink", f"https://drive.google.com/uc?export=download&id={file_id}"),
        "directViewLink": f"https://drive.google.com/uc?export=view&id={file_id}",
        "directDownloadLink": f"https://drive.google.com/uc?export=download&id={file_id}",
        "localPath": str(local_path),
    }


def upload_outputs_to_drive(
    credentials_path: Path, token_path: Path, folder_id: str, upload_paths: list[Path]
) -> list[dict[str, str]]:
    service = build_drive_service(credentials_path=credentials_path, token_path=token_path)
    results: list[dict[str, str]] = []
    for path in upload_paths:
        if path.exists():
            results.append(upload_file_to_drive(service=service, local_path=path, folder_id=folder_id))
    return results


def ensure_drive_subfolder(
    service: object, parent_folder_id: str, subfolder_name: str
) -> str:
    _, HttpError, _, _, _, _, _ = _drive_imports()
    escaped = _drive_query_escape(subfolder_name)
    query = (
        f"'{parent_folder_id}' in parents and "
        f"name = '{escaped}' and "
        "mimeType = 'application/vnd.google-apps.folder' and trashed = false"
    )
    listed = (
        service.files()
        .list(
            q=query,
            fields="files(id,name)",
            pageSize=1,
            supportsAllDrives=True,
            includeItemsFromAllDrives=True,
        )
        .execute()
    )
    files = listed.get("files", [])
    if files:
        return files[0]["id"]

    try:
        created = (
            service.files()
            .create(
                body={
                    "name": subfolder_name,
                    "mimeType": "application/vnd.google-apps.folder",
                    "parents": [parent_folder_id],
                },
                fields="id,name",
                supportsAllDrives=True,
            )
            .execute()
        )
        return created["id"]
    except HttpError as exc:
        err = str(exc)
        if "storageQuotaExceeded" in err:
            raise RuntimeError(
                "Google Drive 하위 폴더 생성 실패: 서비스 계정 quota 제한입니다. "
                "OAuth credentials 또는 Shared Drive 사용이 필요합니다."
            ) from exc
        raise


def generate_gnuplot_figures(combined_csv: Path, plots_dir: Path, pressures: list[float]) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)

    p1, p2, p3 = pressures
    pressure_terms = [
        (p1, "#1f77b4", "0.1 atm" if abs(p1 - 0.1) < 1e-12 else f"{p1:g} atm"),
        (p2, "#d62728", "1 atm" if abs(p2 - 1.0) < 1e-12 else f"{p2:g} atm"),
        (p3, "#2ca02c", "4 atm" if abs(p3 - 4.0) < 1e-12 else f"{p3:g} atm"),
    ]

    scripts: list[tuple[str, int, str, str, Path]] = [
        ("rho_vs_T.gnuplot", 3, "Density", "rho [kg/m^3]", plots_dir / "rho_vs_T.png"),
        ("h_vs_T.gnuplot", 4, "Specific Enthalpy", "h [J/kg]", plots_dir / "h_vs_T.png"),
        ("cp_vs_T.gnuplot", 5, "Specific Heat at Constant Pressure", "Cp [J/(kg K)]", plots_dir / "cp_vs_T.png"),
        (
            "cp_vs_T_0_30000K.gnuplot",
            5,
            "Specific Heat at Constant Pressure (0-30000 K)",
            "Cp [J/(kg K)]",
            plots_dir / "cp_vs_T_0_30000K.png",
        ),
    ]

    for script_name, col, title, ylabel, png_path in scripts:
        x_range_line = ""
        if "0_30000" in script_name:
            x_range_line = "set xrange [0:30000]\n"

        plot_lines = []
        for p, color, label in pressure_terms:
            plot_lines.append(
                (
                    f"'{combined_csv}' u 1:(abs($2-{p:.16g})<1e-12 ? ${col} : 1/0) "
                    f"w l lw 2 lc rgb '{color}' title '{label}'"
                )
            )

        script_text = (
            "set datafile separator ','\n"
            "set terminal pngcairo size 1400,900 enhanced\n"
            f"set output '{png_path}'\n"
            "set grid\n"
            "set key outside right top\n"
            "set xlabel 'Temperature [K]'\n"
            f"set ylabel '{ylabel}'\n"
            f"set title '{title}'\n"
            f"{x_range_line}"
            "set format y '%.3e'\n"
            "plot \\\n"
            + ", \\\n".join(plot_lines)
            + "\n"
        )

        script_path = plots_dir / script_name
        script_path.write_text(script_text, encoding="utf-8")
        subprocess.run(["gnuplot", str(script_path)], check=True)


def main() -> None:
    args = parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / "plots"

    temperature_pf, q_by_species, eion_by_species_eV = load_partition_data(
        args.partition_json.resolve()
    )
    temperatures_eq, equilibrium_by_pressure, pressures = load_equilibrium_data(
        args.equilibrium_csv.resolve()
    )

    if len(pressures) != 3:
        raise ValueError(
            f"Expected 3 pressure cases, found {len(pressures)}: {pressures}"
        )

    masses = {
        "Ar": MASS_AR,
        "Ar+": MASS_AR - MASS_E,
        "Ar2+": MASS_AR - 2.0 * MASS_E,
        "Ar3+": MASS_AR - 3.0 * MASS_E,
        "Ar4+": MASS_AR - 4.0 * MASS_E,
        "e-": MASS_E,
    }

    h_species_j_pf, dlnq_dt_pf = compute_species_enthalpy_j_per_particle(
        temperature_pf=temperature_pf,
        q_by_species=q_by_species,
        eion_by_species_eV=eion_by_species_eV,
    )

    combined_rows: list[dict[str, float]] = []
    per_pressure_rows: dict[float, list[dict[str, float]]] = {}

    for p_atm in pressures:
        state_rows = equilibrium_by_pressure[p_atm]
        temps = [row["T_K"] for row in state_rows]
        rho_list: list[float] = []
        h_list: list[float] = []

        for state in state_rows:
            t = state["T_K"]
            h_lookup = build_h_lookup_for_temperature(
                temperature_pf=temperature_pf,
                h_species_j=h_species_j_pf,
                temperature=t,
            )

            rho = 0.0
            numerator_h = 0.0
            for sp in SPECIES:
                ni = state[sp]
                rho += ni * masses[sp]
                numerator_h += ni * h_lookup[sp]

            if rho <= 0.0:
                raise ValueError(f"Non-positive density at T={t}, P={p_atm} atm")

            h_mix = numerator_h / rho
            rho_list.append(rho)
            h_list.append(h_mix)

        cp_list = numerical_derivative(temps, h_list)

        rows_for_pressure: list[dict[str, float]] = []
        for i, t in enumerate(temps):
            row = {
                "T_K": t,
                "P_atm": p_atm,
                "rho_kg_m3": rho_list[i],
                "h_J_kg": h_list[i],
                "cp_J_kgK": cp_list[i],
            }
            rows_for_pressure.append(row)
            combined_rows.append(row)

        per_pressure_rows[p_atm] = rows_for_pressure

    combined_rows.sort(key=lambda r: (r["P_atm"], r["T_K"]))

    combined_csv = output_dir / "argon_thermo_0p1_1_4atm.csv"
    write_csv(combined_csv, combined_rows, include_pressure=True)

    for p_atm in pressures:
        tag = pressure_tag(p_atm)
        out_csv = output_dir / f"argon_thermo_{tag}atm.csv"
        write_csv(out_csv, per_pressure_rows[p_atm], include_pressure=False)

    cp_peak_rows: list[dict[str, float]] = []
    for p_atm in pressures:
        t_peak, cp_peak = cp_peaks_between(
            per_pressure_rows[p_atm], t_min=CP_PEAK_TMIN, t_max=CP_PEAK_TMAX
        )
        cp_peak_rows.append(
            {
                "P_atm": p_atm,
                "cp_peak_T_K_0_30000": t_peak,
                "cp_peak_J_kgK_0_30000": cp_peak,
            }
        )
    cp_peak_csv = output_dir / "argon_cp_peaks_0_30000K.csv"
    write_cp_peak_csv(cp_peak_csv, cp_peak_rows)
    # Backward-compatible alias to avoid breaking existing links/references.
    legacy_cp_peak_csv = output_dir / "argon_cp_peaks_10000_20000K.csv"
    shutil.copyfile(cp_peak_csv, legacy_cp_peak_csv)

    metadata = {
        "inputs": {
            "partition_json": str(args.partition_json.resolve()),
            "equilibrium_csv": str(args.equilibrium_csv.resolve()),
        },
        "species_order": SPECIES,
        "pressures_atm": pressures,
        "temperature_range_K": {
            "min": min(temperatures_eq),
            "max": max(temperatures_eq),
            "n_points_per_pressure": len(temperatures_eq),
        },
        "constants": {
            "kB_eV_per_K": KB_EV_PER_K,
            "eV_to_J": EV_TO_J,
            "masses_kg": masses,
        },
        "formulas": {
            "rho": "rho = sum_i n_i m_i",
            "H_i_eV": "H_i = (5/2) k_B T + k_B T^2 d(ln Q_int,i)/dT + E_ion,i",
            "H_i_J": "H_i[J] = H_i[eV] * 1.602176634e-19",
            "h": "h = (1/rho) * sum_i n_i H_i",
            "Cp": "Cp = (dh/dT)_P (finite difference)",
        },
        "cp_peak_0_30000K": cp_peak_rows,
        "outputs": {
            "combined_csv": str(combined_csv),
            "cp_peak_csv": str(cp_peak_csv),
            "per_pressure_csv": [
                str(output_dir / f"argon_thermo_{pressure_tag(p)}atm.csv")
                for p in pressures
            ],
            "plots_dir": str(plots_dir),
        },
        "dlnQdT_sample": {
            "Ar_at_300K": dlnq_dt_pf["Ar"][0],
            "Ar_at_15000K": interp_linear(temperature_pf, dlnq_dt_pf["Ar"], 15000.0),
            "e_minus_at_15000K": interp_linear(temperature_pf, dlnq_dt_pf["e-"], 15000.0),
        },
    }
    metadata_path = output_dir / "argon_thermo_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    if not args.skip_plots:
        try:
            generate_gnuplot_figures(
                combined_csv=combined_csv, plots_dir=plots_dir, pressures=pressures
            )
            # Backward-compatible aliases for previously shared plot names.
            legacy_plot_map = {
                plots_dir / "cp_vs_T_0_30000K.png": plots_dir / "cp_vs_T_10000_20000K.png",
                plots_dir / "cp_vs_T_0_30000K.gnuplot": plots_dir / "cp_vs_T_10000_20000K.gnuplot",
            }
            for src, dst in legacy_plot_map.items():
                if src.exists():
                    shutil.copyfile(src, dst)
        except FileNotFoundError:
            print("[WARN] gnuplot not found. Skipping plots.")
        except subprocess.CalledProcessError as exc:
            print(f"[WARN] gnuplot failed ({exc}). Skipping plots.")

    gdrive_manifest_path = None
    if args.upload_to_gdrive:
        folder_id = args.gdrive_folder_id or os.getenv("GDRIVE_FOLDER_ID")
        if not folder_id:
            raise ValueError(
                "Google Drive upload enabled but folder ID is missing. "
                "Set --gdrive-folder-id or GDRIVE_FOLDER_ID."
            )

        target_folder_id = folder_id
        if args.gdrive_subfolder_name:
            drive_service = build_drive_service(
                credentials_path=args.gdrive_credentials.resolve(),
                token_path=args.gdrive_token.resolve(),
            )
            target_folder_id = ensure_drive_subfolder(
                service=drive_service,
                parent_folder_id=folder_id,
                subfolder_name=args.gdrive_subfolder_name,
            )

        upload_paths: list[Path] = [
            combined_csv,
            cp_peak_csv,
            *[output_dir / f"argon_thermo_{pressure_tag(p)}atm.csv" for p in pressures],
        ]
        if not args.skip_plots:
            upload_paths.extend(sorted(plots_dir.glob("*.png")))

        uploaded = upload_outputs_to_drive(
            credentials_path=args.gdrive_credentials.resolve(),
            token_path=args.gdrive_token.resolve(),
            folder_id=target_folder_id,
            upload_paths=upload_paths,
        )
        gdrive_manifest_path = output_dir / "argon_gdrive_upload_manifest.json"
        gdrive_manifest_path.write_text(
            json.dumps(
                {
                    "folder_id": folder_id,
                    "target_folder_id": target_folder_id,
                    "subfolder_name": args.gdrive_subfolder_name,
                    "credentials_path": str(args.gdrive_credentials.resolve()),
                    "uploaded_count": len(uploaded),
                    "uploaded_files": uploaded,
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )

    print(f"[OK] Wrote: {combined_csv}")
    for p in pressures:
        print(f"[OK] Wrote: {output_dir / f'argon_thermo_{pressure_tag(p)}atm.csv'}")
    print(f"[OK] Wrote: {cp_peak_csv}")
    print(f"[OK] Wrote: {metadata_path}")
    if not args.skip_plots:
        print(f"[OK] Wrote plots in: {plots_dir}")
    if gdrive_manifest_path is not None:
        print(f"[OK] Wrote: {gdrive_manifest_path}")
    print(f"[INFO] Pressure cases: {pressures}")
    print(
        f"[INFO] Temperature range: {min(temperatures_eq)} .. {max(temperatures_eq)} K "
        f"({len(temperatures_eq)} points per pressure)"
    )


if __name__ == "__main__":
    main()
