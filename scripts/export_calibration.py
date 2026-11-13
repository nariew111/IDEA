#!/usr/bin/env python3
"""Export a checked-in calibration profile in the format consumed by IDEA-design."""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import json
from pathlib import Path
import shutil


def export_profile(resources: Path, species: str, population: str, density: int,
                   output: Path, historical: bool = False) -> dict:
    manifest = json.loads((resources / "manifest.json").read_text())
    features = set(manifest["fixed_weights"]["retained_features"])
    version = manifest["version"]
    source = (resources / "historical/reader_priorities.tsv" if historical
              else resources / version / "calibration.tsv")
    with source.open() as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t")
                if row["species"] == species and row["population"] == population
                and int(row["density"]) == density]
    if {row["feature"] for row in rows} != features:
        raise ValueError(f"No complete calibration profile for {species}/{population}/{density}")
    grouped = {feature: [row for row in rows if row["feature"] == feature]
               for feature in sorted(features)}
    for feature, records in grouped.items():
        ordered = sorted(records, key=lambda row: int(row["bin"]))
        expected = list(range(1, (2 if feature == "Conserved" else 6) + 1))
        if [int(row["bin"]) for row in ordered] != expected:
            raise ValueError(f"Incomplete or duplicated bins: {feature}")
        column = "historical_value" if historical else "mean_f1"
        values = [Decimal(row[column]) for row in ordered]
        if any(not value.is_finite() or not 0 <= value <= 1 for value in values):
            raise ValueError(f"Invalid calibration value: {feature}")
        if not historical and any(int(row["n_samples"]) != 5 for row in ordered):
            raise ValueError("Reconstruction requires exactly five samples per feature bin")
    output.mkdir(parents=True, exist_ok=True)
    for feature, records in grouped.items():
        directory = output / feature
        directory.mkdir(exist_ok=True)
        ordered = sorted(records, key=lambda row: int(row["bin"]))
        with (directory / "summary.txt").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["INV_ID", *(f"bin_{row['bin']}" for row in ordered)])
            column = "historical_value" if historical else "mean_f1"
            # Do not call a historical last-row fallback a population average.
            label = "historical_effective" if historical else "average"
            writer.writerow([label, *(row[column] for row in ordered)])
    shutil.copyfile(resources / "weights/retained_weights.csv", output / "feature_importance1.csv")
    metadata = {
        "species": species, "population": population, "density": density,
        "calibration": "historical_reader_archive" if historical else version,
        "source_table": source.relative_to(resources).as_posix(),
        "fixed_weights": "unchanged retained nonzero weights from the archived publication object",
        "status": "Historical reader behavior; original ABS not verified" if historical
                  else "New five-sample reconstruction; not the original study ABS",
        "new_biological_performance_evaluated": False,
    }
    (output / "calibration_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, choices=["pig", "human", "chicken"])
    parser.add_argument("--population", required=True)
    parser.add_argument("--density", required=True, type=int,
                        help="Archived representative-chromosome density label, not genome-wide total")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--historical", action="store_true",
                        help="Export actual historical reader priorities instead of new five-sample means")
    parser.add_argument("--resources", type=Path,
                        default=Path(__file__).resolve().parents[1] / "resources/idea_design")
    args = parser.parse_args()
    try:
        metadata = export_profile(args.resources, args.species, args.population, args.density,
                                  args.output_dir, args.historical)
    except (ValueError, FileNotFoundError) as exc:
        parser.error(str(exc))
    print(json.dumps(metadata, indent=2))
    print(f"Use --feature-weight {args.output_dir / 'feature_importance1.csv'}")
    print(f"Use --summary-root {args.output_dir}")


if __name__ == "__main__":
    main()
