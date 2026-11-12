from pathlib import Path

import pandas as pd

from idea.cli import main
from idea.legacy.uldsuite_acc import (
    _load_feature_table_for_selection,
    _select_snps_stratified_by_windows,
)


def test_toy_select(tmp_path):
    root = Path(__file__).resolve().parents[1]
    out = tmp_path / "toy_out"
    main(
        [
            "select",
            "--species",
            "pig",
            "--breeds",
            "toy",
            "--num-bins",
            "3",
            "--data-root",
            str(root / "examples" / "toy"),
            "--output-dir",
            str(out),
        ]
    )
    selected = (out / "toy_selected_snps.txt").read_text().strip().splitlines()
    expected = (root / "examples" / "toy" / "expected_selected_snps.txt").read_text().strip().splitlines()
    assert selected == expected


def test_select_only_accepts_historical_snp_encoded_chromosome():
    root = Path(__file__).resolve().parents[1]
    feature_file = root / "examples" / "toy" / "pig" / "toy" / "toy_feature_8.txt"
    table = _load_feature_table_for_selection(feature_file)
    assert table["CHR"].tolist() == ["1"] * 6


def test_sparse_panel_is_spread_across_all_windows(tmp_path):
    table = pd.DataFrame(
        {
            "SNP": [f"1_{position}" for position in range(10, 101, 10)],
            "POS": list(range(10, 101, 10)),
        }
    )
    selected = _select_snps_stratified_by_windows(
        chrom_df=table,
        weight_dict={},
        summary_root=tmp_path,
        default_feature_bins=6,
        window_bp=10,
        target_snp_count=3,
    )
    assert selected == ["1_10", "1_50", "1_100"]
