import csv
from collections import defaultdict
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from idea.legacy.uldsuite_acc import load_feature_bin_priorities


ROOT = Path(__file__).resolve().parents[1]
RESOURCES = ROOT / "resources/idea_design"


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_archived_weights_and_resource_checksums():
    original = RESOURCES / "weights/feature_importance1.csv"
    assert hashlib.sha256(original.read_bytes()).hexdigest() == (
        "584bfce3791adcfa4384be6d2529635bb9a40ce4df1acf0142a34591f242fda8"
    )
    with (RESOURCES / "weights/retained_weights.csv").open() as handle:
        retained = {row["feature"]: Decimal(row["importance"]) for row in csv.DictReader(handle)}
    assert retained == {"MAF": 1174, "L2": 655, "Conserved": 470, "Conserved_score": 400}
    for record in read_tsv(RESOURCES / "checksums.tsv"):
        assert hashlib.sha256((RESOURCES / record["file"]).read_bytes()).hexdigest() == record["sha256"]


def test_resource_coverage_and_new_calibration_provenance():
    manifest = json.loads((RESOURCES / "manifest.json").read_text())
    assert manifest["coverage"] == {
        "populations": 20, "profiles": 140, "source_summary_files": 560, "feature_bin_rows": 2800
    }
    assert manifest["historical_archive"]["original_five_sample_ABS_verified"] is False
    assert manifest["validation"]["original_manuscript_panels_reproduced"] is False
    membership = defaultdict(set)
    for row in read_tsv(RESOURCES / "reconstructed_20261113/selected_samples.tsv"):
        membership[row["species"], row["population"]].add(row["sample_id"])
        assert row["role"] == "post_study_reconstruction_only"
    assert len(membership) == 20
    assert all(len(ids) == 5 for ids in membership.values())
    groups = defaultdict(list)
    for row in read_tsv(RESOURCES / "reconstructed_20261113/calibration.tsv"):
        groups[row["species"], row["population"], row["density"], row["feature"]].append(row)
        assert int(row["n_samples"]) == 5
        assert Decimal(0) <= Decimal(row["mean_f1"]) <= Decimal(1)
    assert len(groups) == 560
    for key, rows in groups.items():
        expected_count = 2 if key[3] == "Conserved" else 6
        assert sorted(int(row["bin"]) for row in rows) == list(range(1, expected_count + 1))
        assert sorted(int(row["priority"]) for row in rows) == list(range(1, expected_count + 1))


@pytest.mark.parametrize("species,population,density,historical", [
    ("pig", "D", 1000, False),
    ("human", "EUR", 1600, False),
    ("chicken", "RJF", 300, False),
    ("human", "EUR", 1600, True),
])
def test_exported_profile_is_consumed_without_rank_changes(tmp_path, species, population, density, historical):
    command = [sys.executable, str(ROOT / "scripts/export_calibration.py"),
               "--species", species, "--population", population, "--density", str(density),
               "--output-dir", str(tmp_path)]
    if historical:
        command.append("--historical")
    subprocess.run(command, check=True, capture_output=True, text=True)
    table = (RESOURCES / "historical/reader_priorities.tsv" if historical
             else RESOURCES / "reconstructed_20261113/calibration.tsv")
    records = [row for row in read_tsv(table) if row["species"] == species
               and row["population"] == population and int(row["density"]) == density]
    expected = defaultdict(dict)
    for row in records:
        expected[row["feature"]][int(row["bin"])] = int(row["priority"])
    for feature, priority in expected.items():
        actual, count = load_feature_bin_priorities(feature, tmp_path)
        assert actual == priority
        assert count == len(priority)
    metadata = json.loads((tmp_path / "calibration_metadata.json").read_text())
    assert metadata["calibration"] == ("historical_reader_archive" if historical else "reconstructed_20261113")
    assert metadata["new_biological_performance_evaluated"] is False


def test_unknown_profile_is_rejected_without_fallback(tmp_path):
    result = subprocess.run([
        sys.executable, str(ROOT / "scripts/export_calibration.py"),
        "--species", "human", "--population", "EUR", "--density", "999", "--output-dir", str(tmp_path),
    ], text=True, capture_output=True)
    assert result.returncode != 0
    assert "No complete calibration profile" in result.stderr
    assert not list(tmp_path.iterdir())
