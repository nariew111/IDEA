#!/usr/bin/env python3
"""Archive fixed design weights and transparently reconstruct calibration profiles.

Reads historical summary tables supplied by the author. It does not impute
genotypes, retrain weights or reproduce the study's original ABS membership.
Only aggregate bin values and public-reference sample metadata are exported.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import shutil


POPULATIONS = {
    "pig": ["D", "DJ", "DT", "LWU", "MS", "QP", "TC", "TP", "WS", "Y"],
    "human": ["AFR", "AMR", "EAS", "EUR", "SAS"],
    "chicken": ["Glh", "Leg", "RJF", "Tibe", "YpD"],
}
DENSITIES = {
    "pig": [100, 200, 400, 700, 1000, 1500],
    "human": [220, 400, 600, 800, 1200, 1600],
    "chicken": [10, 20, 30, 40, 47, 60, 80, 100, 200, 300],
}
BUILDS = {"pig": "Sscrofa11.1", "human": "GRCh37/hg19", "chicken": "GRCg6a/galGal6"}
REFERENCES = {"pig": "PHARP 4.0", "human": "1000 Genomes / IGSR", "chicken": "ChickenGTEx"}
CHROMOSOMES = {"pig": "6", "human": "22", "chicken": "17"}
WEIGHT_SOURCE = "3_snp_group/results/Species_D/feature_importance1.csv"
WEIGHT_SHA256 = "584bfce3791adcfa4384be6d2529635bb9a40ce4df1acf0142a34591f242fda8"
SEED = "20261113"
SAMPLE_COUNT = 5
VERSION = "reconstructed_20261113"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def choose_samples(sample_ids: set[str], species: str, population: str) -> list[str]:
    """Choose five IDs without consulting imputation outcomes; stable across platforms."""
    if len(sample_ids) < SAMPLE_COUNT:
        raise ValueError(f"Fewer than five common samples: {species}/{population}")

    def key(sample_id: str) -> tuple[str, str]:
        value = f"{SEED}\t{species}\t{population}\t{sample_id}".encode()
        return hashlib.sha256(value).hexdigest(), sample_id

    return sorted(sample_ids, key=key)[:SAMPLE_COUNT]


def write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_summary(path: Path) -> tuple[dict[str, dict], list[str], dict | None, dict]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    if not rows or not fields:
        raise ValueError(f"Empty summary: {path}")
    sample_column = fields[0]
    bins = sorted(
        (field for field in fields if field.startswith("bin_") and field[4:].isdigit()),
        key=lambda field: int(field[4:]),
    )
    if [int(field[4:]) for field in bins] not in [list(range(1, 3)), list(range(1, 7))]:
        raise ValueError(f"Unexpected feature bins: {path}")
    samples = {}
    averages = []
    for row in rows:
        sample_id = row[sample_column].strip()
        if sample_id.lower() == "average":
            averages.append(row)
            continue
        if sample_id in samples:
            raise ValueError(f"Duplicate sample: {path}")
        for column in bins:
            value = Decimal(row[column])
            if not value.is_finite() or not Decimal(0) <= value <= Decimal(1):
                raise ValueError(f"Invalid retained-feature F1: {path}, {column}")
        samples[sample_id] = row
    if len(averages) > 1:
        raise ValueError(f"Multiple Average rows: {path}")
    return samples, bins, averages[0] if averages else None, rows[-1]


def ranks(values: dict[int, Decimal]) -> dict[int, int]:
    ordered = sorted(values, key=lambda index: (-values[index], index))
    return {index: len(values) - position for position, index in enumerate(ordered)}


def build(source_root: Path, output: Path) -> dict:
    source_weights = source_root / WEIGHT_SOURCE
    if digest(source_weights) != WEIGHT_SHA256:
        raise ValueError("Fixed weight source differs from the documented publication object")
    with source_weights.open() as handle:
        weights = list(csv.DictReader(handle))
    retained = [row for row in weights if Decimal(row["importance"]) > 0]
    features = [row["feature"] for row in retained]
    weight_dir = output / "weights"
    weight_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_weights, weight_dir / "feature_importance1.csv")
    with (weight_dir / "retained_weights.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["feature", "importance"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(retained)

    calibration = []
    legacy = []
    sample_rows = []
    provenance = []
    profile_rows = []
    population_rows = []
    for species, populations in POPULATIONS.items():
        for population in populations:
            loaded = {}
            common = None
            for density in DENSITIES[species]:
                for feature in features:
                    relative = Path("3_snp_group") / species / f"{population}_new_{density}" / feature / "summary.txt"
                    parsed = load_summary(source_root / relative)
                    loaded[density, feature] = (relative, parsed)
                    ids = set(parsed[0])
                    common = ids if common is None else common & ids
            assert common is not None
            selected = choose_samples(common, species, population)
            population_rows.append({
                "species": species, "population": population, "genome_build": BUILDS[species],
                "reference_resource": REFERENCES[species], "common_source_samples": len(common),
                "common_sample_ids_sha256": hashlib.sha256("\n".join(sorted(common)).encode()).hexdigest(),
                "selected_samples": selected,
            })
            for position, sample_id in enumerate(selected, 1):
                sample_rows.append({
                    "species": species, "population": population, "selection_order": position,
                    "sample_id": sample_id, "reference_resource": REFERENCES[species],
                    "role": "post_study_reconstruction_only",
                })
            for density in DENSITIES[species]:
                profile_rows.append({
                    "species": species, "population": population, "density": density,
                    "representative_chromosome": CHROMOSOMES[species], "genome_build": BUILDS[species],
                    "n_selected_samples": SAMPLE_COUNT, "features": ",".join(features),
                    "calibration_version": VERSION,
                })
                for feature in features:
                    relative, (samples, bins, average, last) = loaded[density, feature]
                    values = {
                        int(column[4:]): sum(Decimal(samples[sample][column]) for sample in selected) / SAMPLE_COUNT
                        for column in bins
                    }
                    priority = ranks(values)
                    for bin_number, value in values.items():
                        calibration.append({
                            "species": species, "population": population, "density": density,
                            "feature": feature, "bin": bin_number, "mean_f1": str(value),
                            "priority": priority[bin_number], "n_samples": SAMPLE_COUNT,
                        })
                    selected_row = average if average is not None else last
                    row_policy = "explicit_average_row" if average is not None else "last_row_fallback"
                    historical_values = {int(column[4:]): Decimal(selected_row[column]) for column in bins}
                    historical_ranks = ranks(historical_values)
                    for bin_number, value in historical_values.items():
                        legacy.append({
                            "species": species, "population": population, "density": density,
                            "feature": feature, "bin": bin_number, "historical_value": str(value),
                            "priority": historical_ranks[bin_number], "source_row_selection": row_policy,
                        })
                    provenance.append({
                        "species": species, "population": population, "density": density, "feature": feature,
                        "source_file": relative.as_posix(), "source_sha256": digest(source_root / relative),
                        "source_sample_rows": len(samples), "explicit_average_present": average is not None,
                        "historical_reader_row_selection": row_policy,
                    })

    write_tsv(output / VERSION / "calibration.tsv", list(calibration[0]), calibration)
    write_tsv(output / VERSION / "selected_samples.tsv", list(sample_rows[0]), sample_rows)
    write_tsv(output / VERSION / "profiles.tsv", list(profile_rows[0]), profile_rows)
    write_tsv(output / "historical" / "reader_priorities.tsv", list(legacy[0]), legacy)
    write_tsv(output / "source_files.tsv", list(provenance[0]), provenance)
    manifest = {
        "version": VERSION, "created_on": "2026-11-13",
        "status": "post-study reconstruction; not a recovered original ABS",
        "fixed_weights": {
            "source_file": WEIGHT_SOURCE, "sha256": WEIGHT_SHA256,
            "interpretation": "Preserved development-stage empirical scoring object; not a recoverable closed-form SHAP aggregation",
            "retained_features": {row["feature"]: row["importance"] for row in retained},
            "other_feature_weights": "zero; retained_weights.csv omits them without changing scores",
        },
        "sample_selection": {
            "seed": SEED, "count_per_population": SAMPLE_COUNT,
            "scheme": "Ascending SHA256 of UTF-8 seed TAB species TAB population TAB sample_id; ID breaks hash ties",
            "pool": "Intersection of sample IDs in all retained-feature summaries at all listed densities for that population",
            "outcome_used_for_selection": False, "same_samples_at_all_densities": True,
        },
        "aggregation": {
            "method": "Arithmetic mean of five archived sample-level F1 values, separately for each feature bin and density",
            "input_precision": "Archived summaries rounded individual F1 values to three decimal places",
            "priority": "Higher mean F1 receives a larger rank; equal means are ordered by ascending original bin ID",
            "numeric_feature_bin_boundaries": "Not reconstructed or changed; original bin labels are preserved",
        },
        "historical_archive": {
            "interpretation": "Values and ranks actually selected by the historical reader for these source tables",
            "missing_average_behavior": "Historical reader takes the last data row; this is not an estimated population mean",
            "original_five_sample_ABS_verified": False,
        },
        "coverage": {"populations": len(population_rows), "profiles": len(profile_rows),
                     "source_summary_files": len(provenance), "feature_bin_rows": len(calibration)},
        "validation": {"new_imputation_runs": False, "new_prediction_runs": False,
                       "original_manuscript_panels_reproduced": False,
                       "remaining_samples_are_independent_new_validation": False},
        "populations": population_rows,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    generated = sorted(path for path in output.rglob("*") if path.is_file() and path.name not in {"checksums.tsv", "README.md"})
    write_tsv(output / "checksums.tsv", ["file", "sha256"], [
        {"file": path.relative_to(output).as_posix(), "sha256": digest(path)} for path in generated
    ])
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path,
                        help="Author-provided project root containing the historical 3_snp_group tables")
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "resources/idea_design")
    args = parser.parse_args()
    result = build(args.source_root, args.output_dir)
    print(json.dumps(result["coverage"], indent=2))
    print("Saved a dated reconstruction, not the original manuscript ABS.")


if __name__ == "__main__":
    main()
