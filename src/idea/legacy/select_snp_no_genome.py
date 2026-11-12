#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
import re
from pathlib import Path

# ------------------------- 参数解析 -------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Weighted SNP selection using per-feature summary")
    parser.add_argument("--species", required=True, help="Species name, e.g., pig")
    parser.add_argument("--breeds", required=True, help="Space-separated list of breeds, e.g., 'DJ DT D'")
    parser.add_argument("--num-bins", type=int, default=1000, help="Number of position bins for SNP selection")
    parser.add_argument("--data-root", required=True, help="Root directory of prepared feature files")
    parser.add_argument("--output-dir", required=True, help="Directory to save selected SNP files")
    parser.add_argument("--default-feature-bins", type=int, default=6, help="Default bins for feature discretization")
    return parser.parse_args()

# ------------------------- 数据加载 -------------------------
def load_feature_file(feature_file):
    df = pd.read_csv(feature_file, sep="\t", dtype=str)
    if 'SNP' not in df.columns:
        raise ValueError(f"Feature file {feature_file} must contain 'SNP' column")
    for col in df.columns:
        if col == 'SNP':
            continue
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def load_feature_weights(weight_file):
    df = pd.read_csv(weight_file)
    val_col = 'importance' if 'importance' in df.columns else ('weight' if 'weight' in df.columns else df.columns[1])
    weight_dict = dict(zip(df['feature'].astype(str), df[val_col].astype(float)))
    return weight_dict

# ------------------------- summary -> bin 优先级 -------------------------
def find_summary_file(feature, breed_dir):
    candidates = [
        os.path.join(breed_dir, feature, "summary.txt"),
        os.path.join(breed_dir, f"{feature}_summary.txt"),
        os.path.join(breed_dir, f"{feature}.summary.txt"),
        os.path.join(breed_dir, f"{feature}.txt")
    ]
    for p in candidates:
        if os.path.isfile(p):
            return p
    return None

def load_feature_bin_priorities(feature, breed_dir, default_n_bins=6):
    summary_file = find_summary_file(feature, breed_dir)
    if summary_file is None:
        # identity mapping
        return {i: i for i in range(1, default_n_bins+1)}, default_n_bins

    try:
        df = pd.read_csv(summary_file, sep=r'\s+', engine='python', dtype=str)
    except:
        return {i: i for i in range(1, default_n_bins+1)}, default_n_bins

    first_col = df.columns[0]
    mask_avg = df[first_col].astype(str).str.strip().str.lower() == 'average'
    avg_row = df[mask_avg].iloc[0] if mask_avg.any() else df.iloc[-1]

    bin_cols = [c for c in df.columns if re.search(r'bin[_]?\d+$', c, flags=re.I)]
    if not bin_cols:
        bin_cols = [c for c in df.columns if re.match(r'^\d+$', c)]
    if not bin_cols:
        return {i: i for i in range(1, default_n_bins+1)}, default_n_bins

    acc_map = {}
    for c in bin_cols:
        m = re.search(r'(\d+)$', c)
        if not m:
            continue
        bin_idx = int(m.group(1))
        try:
            acc_val = float(avg_row[c])
        except:
            acc_val = 0.0
        acc_map[bin_idx] = acc_val

    n_bins = max(acc_map.keys()) if acc_map else default_n_bins
    sorted_bins = sorted(acc_map.items(), key=lambda x: (x[1], -x[0]), reverse=True)
    priority_map = {}
    rank_val = n_bins
    for bin_idx, _ in sorted_bins:
        priority_map[bin_idx] = rank_val
        rank_val -= 1
    for b in range(1, n_bins+1):
        if b not in priority_map:
            priority_map[b] = max(1, rank_val)
            rank_val -= 1
    return priority_map, n_bins

# ------------------------- SNP BIN -------------------------
def process_position_bins(df, num_bins):
    df['POS'] = df['SNP'].str.split('_').str[-1].astype(int)
    bins = np.linspace(df['POS'].min(), df['POS'].max(), num_bins+1)
    df['SNP_BIN'] = pd.cut(df['POS'], bins, labels=False, include_lowest=True) + 1
    mid_snp_map = {}
    for bin_id in range(num_bins):
        start = bins[bin_id]
        end = bins[bin_id+1]
        if bin_id == num_bins - 1:
            sub_df = df[(df['POS'] >= start) & (df['POS'] <= end)].copy()
        else:
            sub_df = df[(df['POS'] >= start) & (df['POS'] < end)].copy()
        if not sub_df.empty:
            midpoint = (start + end) / 2
            sub_df['dist'] = (sub_df['POS'] - midpoint).abs()
            mid_snp = sub_df.nsmallest(1, 'dist')['SNP'].values[0]
            mid_snp_map[bin_id+1] = mid_snp
    return mid_snp_map

# ------------------------- SNP选择 -------------------------
def select_snps(df, weight_dict, breed_dir, default_feature_bins, mid_snp_map):
    feature_meta = {}
    for feat, w in weight_dict.items():
        priority_map, n_bins = load_feature_bin_priorities(feat, breed_dir, default_n_bins=default_feature_bins)
        feature_meta[feat] = {'weight': float(w), 'priority_map': priority_map, 'n_bins': n_bins}
        feat_bin_col = f"{feat}_BIN"
        if feat_bin_col not in df.columns and feat in df.columns:
            bins = n_bins if feat.lower() != 'conserved' else 2
            df[feat_bin_col] = pd.cut(df[feat].fillna(0), bins=bins, labels=False)+1
        elif feat_bin_col not in df.columns:
            df[feat_bin_col] = 0

    selected_snps = []
    for snp_bin in sorted(df['SNP_BIN'].unique()):
        candidates = df[df['SNP_BIN']==snp_bin].copy()
        if candidates.empty:
            selected_snps.append(mid_snp_map[snp_bin])
            continue
        weighted_score = np.zeros(len(candidates))
        for feat, meta in feature_meta.items():
            w = meta['weight']
            pmap = meta['priority_map']
            bins_array = candidates[f"{feat}_BIN"].fillna(0).astype(int).values
            for bin_k, priority in pmap.items():
                mask = (bins_array == bin_k)
                if mask.any():
                    weighted_score[mask] += priority*w
        candidates['weighted_score'] = weighted_score
        max_score = candidates['weighted_score'].max()
        top_candidates = candidates[candidates['weighted_score']==max_score].copy()
        mid_pos = int(mid_snp_map.get(snp_bin, top_candidates['POS'].median()))
        top_candidates['dist'] = (top_candidates['POS'] - mid_pos).abs()
        best_snp = top_candidates.nsmallest(1, 'dist')['SNP'].values[0]
        selected_snps.append(best_snp)
    return selected_snps

# ------------------------- 保存 -------------------------
def save_snp_list(output_dir, breed, snp_list):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    snp_file = output_dir / f"{breed}_selected_snps.txt"
    pd.Series(snp_list).to_csv(snp_file, index=False, header=False)
    print(f"[INFO] Saved {breed} SNPs to {snp_file}")
    return snp_file

# ------------------------- 主流程 -------------------------
def main():
    args = parse_args()
    breeds = args.breeds.split()
    for breed in breeds:
        print(f"[INFO] Processing breed {breed}")
        breed_dir = os.path.join(args.data_root, args.species, breed)
        feature_file = os.path.join(breed_dir, f"{breed}_feature_8.txt")
        weight_file = os.path.join(breed_dir, "feature_importance1.csv")
        df = load_feature_file(feature_file)
        weight_dict = load_feature_weights(weight_file)
        mid_snp_map = process_position_bins(df, args.num_bins)
        selected_snps = select_snps(df, weight_dict, breed_dir, args.default_feature_bins, mid_snp_map)
        save_snp_list(args.output_dir, breed, selected_snps)

if __name__ == "__main__":
    main()
