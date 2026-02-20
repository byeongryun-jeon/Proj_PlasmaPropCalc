#!/usr/bin/env python3
"""Extract per-database e-Ar momentum-transfer tables from LXCat text export.

Input:
  - LXCat text export file (default: /home/brjeon/matfLookUpPython/lxcatData.txt)

Output:
  - One CSV per database under data/raw/murphy/collision_integrals/lxcat/
  - Summary CSV/JSON with point counts and energy range.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


M2_TO_A2 = 1.0e20


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Extract e-Ar ELASTIC blocks from LXCat text export into per-database CSV files."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("/home/brjeon/matfLookUpPython/lxcatData.txt"),
        help="Path to LXCat text export file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "data" / "raw" / "murphy" / "collision_integrals" / "lxcat",
        help="Directory where per-database CSV files are written.",
    )
    parser.add_argument(
        "--process",
        choices=["ELASTIC", "EFFECTIVE"],
        default="ELASTIC",
        help="LXCat process block type to extract.",
    )
    return parser.parse_args()


def slugify(name: str) -> str:
    s = name.lower()
    s = s.replace("&", "and")
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def parse_two_floats(line: str) -> tuple[float, float] | None:
    parts = [p for p in re.split(r"[,\s\t]+", line.strip()) if p]
    nums: list[float] = []
    for p in parts:
        try:
            nums.append(float(p))
        except ValueError:
            continue
    if len(nums) < 2:
        return None
    return nums[0], nums[1]


def parse_blocks(text: str, process: str) -> list[dict[str, object]]:
    lines = text.splitlines()
    blocks: list[dict[str, object]] = []
    current_db = ""
    i = 0
    process_key = process.strip().upper()

    while i < len(lines):
        s = lines[i].strip()

        if s.startswith("DATABASE:"):
            current_db = s.split(":", 1)[1].strip()
            i += 1
            continue

        if s.upper() != process_key:
            i += 1
            continue

        target = lines[i + 1].strip() if i + 1 < len(lines) else ""
        j = i + 2
        species_ok = False
        while j < len(lines) and not re.match(r"^-{5,}$", lines[j].strip()):
            low = lines[j].strip().lower().replace("  ", " ")
            if low.startswith("species:") and ("e / ar" in low or "e/ar" in low):
                species_ok = True
            j += 1
            if j - i > 500:
                break

        if j >= len(lines) or not re.match(r"^-{5,}$", lines[j].strip()):
            i += 1
            continue

        j += 1
        pts: list[tuple[float, float]] = []
        while j < len(lines) and not re.match(r"^-{5,}$", lines[j].strip()):
            parsed = parse_two_floats(lines[j])
            if parsed is not None:
                e, q = parsed
                if e >= 0.0 and q > 0.0:
                    pts.append((e, q))
            j += 1

        if pts and target.startswith("Ar") and species_ok:
            blocks.append(
                {
                    "database": current_db if current_db else "unknown",
                    "target": target,
                    "species_ok": species_ok,
                    "process": process_key,
                    "points": pts,
                }
            )

        i = j + 1

    return blocks


def pick_best_by_database(blocks: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    for b in blocks:
        db = str(b["database"])
        n = len(b["points"])  # type: ignore[arg-type]
        prev = best.get(db)
        if prev is None or n > len(prev["points"]):  # type: ignore[arg-type]
            best[db] = b
    return best


def write_model_csv(path: Path, db: str, process: str, target: str, points: list[tuple[float, float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=["E_eV", "Qm_A2", "source", "database", "process", "target"],
        )
        writer.writeheader()
        src = f"LXCat:{db}:{process}"
        for e, q_m2 in sorted(points, key=lambda p: p[0]):
            writer.writerow(
                {
                    "E_eV": f"{e:.12g}",
                    "Qm_A2": f"{q_m2 * M2_TO_A2:.12g}",
                    "source": src,
                    "database": db,
                    "process": process,
                    "target": target,
                }
            )


def main() -> None:
    args = parse_args()
    in_path = args.input.resolve()
    out_dir = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    text = in_path.read_text(encoding="utf-8", errors="ignore")
    blocks = parse_blocks(text, process=args.process)
    by_db = pick_best_by_database(blocks)

    summary_rows: list[dict[str, object]] = []
    for db in sorted(by_db.keys()):
        block = by_db[db]
        pts = sorted(block["points"], key=lambda p: p[0])  # type: ignore[arg-type]
        slug = slugify(db)
        out_csv = out_dir / f"lxcat_{slug}_e_ar_{args.process.lower()}.csv"
        write_model_csv(
            out_csv,
            db=db,
            process=str(block["process"]),
            target=str(block["target"]),
            points=pts,  # type: ignore[arg-type]
        )
        summary_rows.append(
            {
                "database": db,
                "slug": slug,
                "process": block["process"],
                "target": block["target"],
                "n_points": len(pts),
                "E_min_eV": pts[0][0],
                "E_max_eV": pts[-1][0],
                "output_csv": str(out_csv),
            }
        )
        print(f"[OK] {db}: {len(pts)} points -> {out_csv}")

    summary_csv = out_dir / f"lxcat_e_ar_{args.process.lower()}_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=[
                "database",
                "slug",
                "process",
                "target",
                "n_points",
                "E_min_eV",
                "E_max_eV",
                "output_csv",
            ],
        )
        writer.writeheader()
        for r in summary_rows:
            writer.writerow(r)

    summary_json = out_dir / f"lxcat_e_ar_{args.process.lower()}_summary.json"
    summary_json.write_text(
        json.dumps(
            {
                "input_file": str(in_path),
                "process": args.process,
                "n_databases": len(summary_rows),
                "models": summary_rows,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print(f"[OK] Wrote: {summary_csv}")
    print(f"[OK] Wrote: {summary_json}")


if __name__ == "__main__":
    main()

