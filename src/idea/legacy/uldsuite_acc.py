#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED, as_completed
import glob
import re
import time
import traceback
import yaml
import logging
try:
    from filelock import FileLock
except ModuleNotFoundError:
    class FileLock:
        def __init__(self, *args, **kwargs):
            raise RuntimeError(
                "Missing dependency 'filelock'. Install the release environment with "
                "`conda env create -f environment.yml` or `pip install filelock`."
            )

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    def tqdm(iterable=None, *args, **kwargs):
        return iterable if iterable is not None else []
import shutil

def setup_logger(name, log_file, level=logging.DEBUG): 
    """Function to setup as many loggers as you want""" 
    log_file = Path(log_file) 
    log_file.parent.mkdir(parents=True, exist_ok=True) 

    # Clear existing handlers 
    logger = logging.getLogger(name) 
    if logger.hasHandlers(): 
        logger.handlers.clear() 

    logger.setLevel(level) 

    # File handler 
    fh = logging.FileHandler(log_file, mode='w') 
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')) 
    logger.addHandler(fh) 

    # Console handler 
    ch = logging.StreamHandler(sys.stdout) 
    ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s')) 
    logger.addHandler(ch) 

    return logger 

def run_cmd(cmd_list, logger, is_bcftools=False, bcftools_path=None, cwd=None, extra_env=None): 
    """Runs a command, logs output, and handles errors robustly.""" 
    cmd_str = ' '.join(map(str, cmd_list)) 
    logger.info(f"Running command: {cmd_str}") 
    if cwd:
        logger.debug(f"Working directory: {cwd}")

    env = os.environ.copy() 
    if is_bcftools and bcftools_path: 
        # Resolve full path if it's just a command name or relative path
        if not Path(bcftools_path).is_absolute():
            resolved_path = shutil.which(bcftools_path)
            if resolved_path:
                bcftools_path = resolved_path

        # Add parent of bcftools to LD_LIBRARY_PATH if not already there 
        bcftools_dir = Path(bcftools_path).parent 
        # Also add the ../lib directory which is common for custom builds
        bcftools_lib_dir = bcftools_dir.parent / 'lib'
        
        # Manually check common lib paths relative to conda env
        possible_lib_dirs = [
            bcftools_lib_dir,
            bcftools_dir.parent / 'lib',
            bcftools_dir / '../lib'
        ]

        current_ld = env.get('LD_LIBRARY_PATH', '')
        new_ld_parts = []
        
        for lib_d in possible_lib_dirs:
            resolved = lib_d.resolve()
            if resolved.exists() and str(resolved) not in current_ld.split(':'):
                 new_ld_parts.append(str(resolved))

        if str(bcftools_dir) not in current_ld.split(':'):
            new_ld_parts.append(str(bcftools_dir))
        
        if new_ld_parts:
            env['LD_LIBRARY_PATH'] = f"{':'.join(new_ld_parts)}:{current_ld}"
            logger.debug(f"Updated LD_LIBRARY_PATH: {env['LD_LIBRARY_PATH']}") 
            
        # Update PATH so child processes can find bcftools.
        # Append (not prepend) to avoid overriding the active Python env.
        current_path = env.get('PATH', '')
        if str(bcftools_dir) not in current_path.split(':'):
            env['PATH'] = f"{current_path}:{bcftools_dir}" if current_path else str(bcftools_dir)
            logger.debug(f"Updated PATH: {env['PATH']}") 

    if extra_env:
        for k, v in extra_env.items():
            env[str(k)] = str(v)

    try: 
        # Using subprocess.run with binary capture and manual decoding to handle non-UTF8 output
        result = subprocess.run( 
            cmd_list, 
            capture_output=True, 
            text=False, # Capture raw bytes
            stdin=subprocess.DEVNULL,
            env=env, 
            cwd=cwd,
            check=False 
        ) 
        
        # Manually decode stdout/stderr
        result.stdout = result.stdout.decode('utf-8', errors='replace') if result.stdout else ""
        result.stderr = result.stderr.decode('utf-8', errors='replace') if result.stderr else ""

        if result.stdout: 
            logger.debug(f"STDOUT:\n{result.stdout.strip()}") 
        if result.stderr: 
            # Log stderr as debug, unless return code is non-zero 
            if result.returncode == 0: 
                logger.debug(f"STDERR:\n{result.stderr.strip()}") 
            else: 
                logger.error(f"STDERR:\n{result.stderr.strip()}") 

        # Raise an exception if the command failed 
        if result.returncode != 0: 
            raise subprocess.CalledProcessError( 
                result.returncode, cmd_str, output=result.stdout, stderr=result.stderr 
            ) 
        
        logger.info(f"Command finished successfully: {cmd_str}") 
        return result 

    except FileNotFoundError: 
        logger.error(f"Command not found: {cmd_list[0]}") 
        raise 
    except subprocess.CalledProcessError as e: 
        logger.error(f"Command failed with exit code {e.returncode}: {cmd_str}") 
        # No need to log stderr again as it's already logged above 
        raise 
    except Exception as e: 
        logger.error(f"An unexpected error occurred while running command: {cmd_str}", exc_info=True) 
        raise

def parse_args():
    # First pass: Parse only --config to load defaults, disable help to prevent early exit
    parser = argparse.ArgumentParser(description="Weighted SNP selection and imputation (with per-feature summary-based bin priorities)", add_help=False)
    parser.add_argument("--config", help="Path to config YAML file")
    parser.add_argument("-h", "--help", action="store_true", help="show this help message and exit")

    # Dynamically add arguments from a config file or set as required if no config
    args, unknown = parser.parse_known_args()
    select_only_requested = ("--select-only" in unknown)
    config = {}
    if args.config and os.path.exists(args.config):
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)

    def get_config_val(key, default=None, required=False):
        val = config.get(key, default)
        if required and val is None:
            # If help is requested, we don't want to crash on missing required args yet
            if args.help:
                return None
            raise ValueError(f"Configuration for '{key}' is missing in YAML and has no default.")
        return val

    # Check if config is loaded to determine requirement status
    is_config_loaded = bool(config)

    # Define arguments and set defaults from config
    require_main = (not is_config_loaded) and (not args.help) and (not select_only_requested)

    parser.add_argument("--feature-file", default=get_config_val('feature_file'), required=require_main, help="Path to feature file")
    parser.add_argument("--feature-weight", default=get_config_val('feature_weight'), required=require_main, help="Path to feature importance CSV")
    parser.add_argument("--summary-root", default=get_config_val('summary_root'), required=require_main, help="Root dir for feature summaries")
    parser.add_argument("--reference-vcf", default=get_config_val('reference_vcf'), required=require_main, help="Reference VCF path")
    parser.add_argument("--output-dir", default=get_config_val('output_dir', './imputation_output'), help="Base output directory")
    parser.add_argument("--beagle-jar", default=get_config_val('beagle_jar'), required=require_main, help="Path to Beagle JAR")
    parser.add_argument("--conform-jar", default=get_config_val('conform_jar'), required=require_main, help="Path to conform-gt JAR")
    parser.add_argument("--accuracy-script", default=get_config_val('accuracy_script'), required=require_main, help="Accuracy calculation script")
    parser.add_argument("--benchmark", default=get_config_val('benchmark'), required=require_main, help="Path to imputation benchmark file")

    parser.add_argument("--plink", default=get_config_val('plink'), required=require_main, help="Path to plink2 executable")
    parser.add_argument("--bcftools", default=get_config_val('bcftools'), required=require_main, help="Path to bcftools executable")
    parser.add_argument("--chr", type=str, default=get_config_val('chr', None), required=False, help="Chromosome(s) to process. Examples: '17' or '1,2,3' or '1-18'. If omitted, process all chromosomes present in the feature file.")
    parser.add_argument("--test-chr", type=str, default=get_config_val('test_chr', None), required=False, help="Chromosome to run the imputation/accuracy test on (e.g. '17'). Special value: 'shortest' selects the chromosome with the smallest max POS from the feature file. If omitted, test all chromosomes in --chr (or all in feature file).")

    parser.add_argument("--threads", type=int, default=get_config_val('threads', 4), help="Number of parallel threads")
    parser.add_argument("--beagle-mem", type=int, default=get_config_val('beagle_mem', 200), help="Memory for Beagle in GB")
    parser.add_argument("--num-bins", type=int, default=get_config_val('num_bins', 1000), help="Number of SNP position bins")
    parser.add_argument("--snp-bin-size", type=int, default=get_config_val('snp_bin_size', 50000), help="Size of SNP bins in base pairs for global selection.")
    parser.add_argument("--target-snp-count", type=int, default=get_config_val('target_snp_count', 1000), help="Target number of SNPs to select per chromosome.")
    parser.add_argument("--target-snp-total", type=int, default=get_config_val('target_snp_total', None), help="Target number of SNPs to select across all chromosomes (overrides --target-snp-count).")
    parser.add_argument("--default-feature-bins", type=int, default=get_config_val('default_feature_bins', 6), help="Default number of feature bins")
    parser.add_argument("--max-retries", type=int, default=get_config_val('max_retries', 2), help="Max retries for failed tasks")
    parser.add_argument("--window-size", type=int, default=get_config_val('window_size', 5000000), help="Genomic window size for chunking")
    parser.add_argument("--keep-intermediate", action='store_true', default=get_config_val('keep_intermediate', False), help="Keep intermediate files")
    parser.add_argument("--cleanup", action='store_true', default=get_config_val('cleanup', True), help="Clean up chunk directories after processing")
    parser.add_argument("--skip-imputation", action='store_true', default=get_config_val('skip_imputation', False), help="Skip the imputation step and merge phased VCFs directly.")

    # Dynamic window arguments
    parser.add_argument("--dynamic-window", action='store_true', default=get_config_val('dynamic_window', False), help="Enable dynamic window sizing based on SNP density and count.")
    parser.add_argument("--max-snps-per-chunk", type=int, default=get_config_val('max_snps_per_chunk', 15000), help="Maximum number of SNPs per chunk.")
    parser.add_argument("--max-chunk-size", type=int, default=get_config_val('max_chunk_size', 1000000), help="Maximum physical size of a genomic chunk in base pairs.")

    # New argument for enabling chunking (default is NO chunking)
    parser.add_argument("--enable-chunk", action='store_true', default=get_config_val('enable_chunk', False), help="Enable chunking logic. Default is False (process whole chromosome).")

    # Resource management arguments
    parser.add_argument("--max-total-mem", type=int, default=get_config_val('max_total_mem', 200), help="Maximum total memory (GB) for all concurrent Beagle tasks.")
    parser.add_argument("--max-concurrent-tasks", type=int, default=get_config_val('max_concurrent_tasks', 5), help="Maximum number of concurrent imputation tasks.")
    parser.add_argument("--select-only", action="store_true", default=get_config_val('select_only', False), help="Only perform SNP selection and write a SNP list, then exit.")
    parser.add_argument("--select-out", type=str, default=get_config_val('select_out', None), help="Output path for selected SNP list (one SNP ID per line).")

    if args.help:
        parser.print_help()
        sys.exit(0)

    # Re-parse arguments now that they are all defined
    return parser.parse_args()

def _load_feature_table_for_selection(feature_file):
    df = pd.read_csv(feature_file, sep="\t")
    if 'SNP' not in df.columns:
        raise ValueError("feature file must contain 'SNP' column")
    if 'CHR' not in df.columns:
        # Historical IDEA feature tables encode chromosome and position in the
        # SNP identifier (for example ``1_12345``) and do not carry a separate
        # CHR column. Keep that release format valid for the select-only path.
        parsed_chr = df['SNP'].astype(str).str.replace('^chr', '', regex=True).str.split('_').str[0]
        if parsed_chr.eq('').any():
            raise ValueError("feature file must contain 'CHR' column or chromosome-prefixed SNP IDs")
        df['CHR'] = parsed_chr

    df['SNP'] = df['SNP'].astype(str)
    df['CHR'] = df['CHR'].astype(str).str.replace('^chr', '', regex=True)

    if 'POS' not in df.columns:
        df['POS'] = df['SNP'].apply(lambda x: int(str(x).split('_')[-1]))
    else:
        df['POS'] = pd.to_numeric(df['POS'], errors='coerce')

    for col in df.columns:
        if col in ('SNP', 'CHR'):
            continue
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df.dropna(subset=['POS'])
    df['POS'] = df['POS'].astype(int)
    return df

def _select_snps_stratified_by_windows(chrom_df, weight_dict, summary_root, default_feature_bins, window_bp, target_snp_count):
    if chrom_df.empty:
        return []
    df = chrom_df.copy()

    feature_meta = {}
    for feat, w in weight_dict.items():
        priority_map, n_bins = load_feature_bin_priorities(feat, summary_root, default_n_bins=default_feature_bins)
        feature_meta[feat] = {'weight': float(w), 'priority_map': priority_map, 'n_bins': n_bins}
        feat_bin_col = f"{feat}_BIN"
        if feat_bin_col not in df.columns:
            if feat in df.columns:
                target_bins = n_bins if feat.lower() != 'conserved' else 2
                series_data = df[feat].fillna(0)
                if feat.lower() == 'conserved':
                    df[feat_bin_col] = pd.cut(series_data, bins=target_bins, labels=False) + 1
                else:
                    try:
                        qcut_res = pd.qcut(series_data, q=target_bins, labels=False, duplicates='drop')
                        if qcut_res.nunique() < target_bins:
                            df[feat_bin_col] = pd.cut(series_data, bins=target_bins, labels=False) + 1
                        else:
                            df[feat_bin_col] = qcut_res + 1
                    except ValueError:
                        df[feat_bin_col] = pd.cut(series_data, bins=target_bins, labels=False) + 1
            else:
                df[feat_bin_col] = 0

    df['weighted_score'] = np.zeros(len(df))
    for feat, meta in feature_meta.items():
        w = meta['weight']
        pmap = meta['priority_map']
        feat_bin_col = f"{feat}_BIN"
        if feat_bin_col in df.columns:
            for bin_k, priority in pmap.items():
                if priority == 0:
                    continue
                mask = (df[feat_bin_col].fillna(0).astype(int) == int(bin_k))
                df.loc[mask, 'weighted_score'] += priority * w

    df['WIN'] = (df['POS'] // int(window_bp)).astype(int)
    windows = df['WIN'].drop_duplicates().sort_values().tolist()
    if not windows:
        return []

    target_snp_count = min(int(target_snp_count), len(df))
    per_win = target_snp_count // len(windows)
    remainder = target_snp_count % len(windows)
    extra_windows = set()
    if remainder:
        # Spread the remainder across the full chromosome. The old
        # max(1, target // windows) implementation selected one SNP from every
        # window and then truncated the list, concentrating sparse panels at
        # the start whenever target < number of windows.
        extra_windows = {
            int(index)
            for index in np.linspace(0, len(windows) - 1, num=remainder, dtype=int)
        }
    selected = []

    for window_index, w_id in enumerate(windows):
        sub = df[df['WIN'] == w_id].sort_values(by=['weighted_score', 'POS'], ascending=[False, True])
        take = min(per_win + (1 if window_index in extra_windows else 0), len(sub))
        if take > 0:
            selected.extend(sub['SNP'].head(take).tolist())

    if len(selected) < int(target_snp_count):
        pool = df[~df['SNP'].isin(selected)].sort_values(by=['weighted_score', 'POS'], ascending=[False, True])
        need = int(target_snp_count) - len(selected)
        selected.extend(pool['SNP'].head(need).tolist())

    return selected[:int(target_snp_count)]

def _allocate_target_counts(chrom_to_candidate_count, total_target):
    chrom_items = [(int(k), int(v)) for k, v in chrom_to_candidate_count.items() if int(v) > 0]
    if total_target is None:
        return {}
    total_target = int(total_target)
    if total_target <= 0:
        return {chrom: 0 for chrom, _ in chrom_items}
    if not chrom_items:
        return {}

    denom = sum(v for _, v in chrom_items)
    if denom <= 0:
        return {chrom: 0 for chrom, _ in chrom_items}

    raw = []
    for chrom, n in chrom_items:
        r = (total_target * n) / denom
        base = int(np.floor(r))
        raw.append((chrom, n, r, base, r - base))

    alloc = {chrom: min(base, n) for chrom, n, _, base, _ in raw}
    used = sum(alloc.values())
    remaining = total_target - used

    if remaining <= 0:
        return alloc

    raw_sorted = sorted(raw, key=lambda x: (x[4], x[2], x[1], -x[0]), reverse=True)
    for chrom, n, _, _, _ in raw_sorted:
        if remaining <= 0:
            break
        if alloc.get(chrom, 0) < n:
            alloc[chrom] = alloc.get(chrom, 0) + 1
            remaining -= 1

    return alloc

def _parse_chr_arg(chr_arg, feature_df=None):
    if chr_arg is None or str(chr_arg).strip() == "":
        if feature_df is None or 'CHR' not in feature_df.columns:
            return []
        chrom_series = pd.to_numeric(feature_df['CHR'], errors='coerce').dropna()
        return sorted({int(x) for x in chrom_series.tolist() if int(x) > 0})

    chr_arg = str(chr_arg).strip()
    if '-' in chr_arg:
        start_s, end_s = chr_arg.split('-', 1)
        return list(range(int(start_s), int(end_s) + 1))
    if ',' in chr_arg:
        return [int(c) for c in chr_arg.split(',') if c.strip() != ""]
    return [int(chr_arg)]

def _choose_test_chromosome(test_chr_arg, feature_df, selection_chromosomes):
    if test_chr_arg is None or str(test_chr_arg).strip() == "":
        return None

    test_chr_arg = str(test_chr_arg).strip().lower()
    if test_chr_arg == "shortest":
        if feature_df is None or feature_df.empty or 'CHR' not in feature_df.columns or 'POS' not in feature_df.columns:
            return None
        lengths = []
        for chrom in selection_chromosomes:
            sub = feature_df[feature_df['CHR'] == chrom]
            if sub.empty:
                continue
            max_pos = int(pd.to_numeric(sub['POS'], errors='coerce').dropna().max()) if not sub['POS'].empty else 0
            if max_pos > 0:
                lengths.append((max_pos, int(chrom)))
        if not lengths:
            return None
        lengths.sort()
        return lengths[0][1]

    return int(test_chr_arg)

def run_select_only(args):
    feature_df = _load_feature_table_for_selection(args.feature_file)
    weight_dict = load_feature_weights(args.feature_weight)

    chr_list = _parse_chr_arg(args.chr, feature_df=feature_df)

    chrom_to_n = {}
    for chrom in chr_list:
        chrom_str = str(chrom)
        chrom_to_n[chrom] = int((feature_df['CHR'] == chrom_str).sum())
    alloc = _allocate_target_counts(chrom_to_n, getattr(args, "target_snp_total", None))
    selected_all = []
    for chrom in chr_list:
        chrom_str = str(chrom)
        chrom_df = feature_df[feature_df['CHR'] == chrom_str]
        per_chrom_target = int(alloc.get(chrom, args.target_snp_count)) if alloc else int(args.target_snp_count)
        snps = _select_snps_stratified_by_windows(
            chrom_df=chrom_df,
            weight_dict=weight_dict,
            summary_root=args.summary_root,
            default_feature_bins=args.default_feature_bins,
            window_bp=args.snp_bin_size,
            target_snp_count=per_chrom_target
        )
        selected_all.extend(snps)

    # A panel is a set, but its serialized representation must be canonical for
    # reproducible manifests. Preserve the feature table's genomic order and
    # drop any accidental duplicates from chromosome allocation.
    selected_set = set(map(str, selected_all))
    selected_all = [snp for snp in feature_df['SNP'].astype(str) if snp in selected_set]

    out_path = args.select_out
    if out_path is None:
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = str(out_dir / "selected_snps.txt")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    pd.Series(selected_all).to_csv(out_path, index=False, header=False)
    print(f"Saved selected SNPs: {out_path}")

# ------------------------- 数据加载 -------------------------
def load_features(feature_file):
    df = pd.read_csv(feature_file, sep="\t")
    # 尝试把 feature 列转换成 numeric（除了 SNP 列）
    if 'SNP' not in df.columns:
        raise ValueError("feature file must contain 'SNP' column")

    # 标准化 SNP 列中的染色体名称（移除 'chr' 前缀）
    df['SNP'] = df['SNP'].astype(str).str.replace('^chr', '', regex=True)
    # Add POS column by parsing SNP column
    df['POS'] = df['SNP'].apply(lambda x: int(x.split('_')[-1]))

    for col in df.columns:
        if col == 'SNP':
            continue
        # 尝试转换为 numeric，非数值设为 NaN
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def load_feature_weights(weight_file):
    df = pd.read_csv(weight_file)
    # 支持列名 'importance' 或 'weight'
    if 'importance' in df.columns:
        val_col = 'importance'
    elif 'weight' in df.columns:
        val_col = 'weight'
    else:
        # 如果只有两列，第二列作为权重
        if df.shape[1] >= 2:
            val_col = df.columns[1]
        else:
            raise ValueError("feature weight file must contain feature and importance/weight column")
    weight_dict = dict(zip(df['feature'].astype(str), df[val_col].astype(float)))
    return weight_dict

# ------------------------- per-feature summary -> bin priority -------------------------
def find_summary_file_for_feature(feature, summary_root):
    """尝试在 summary_root 下为 feature 找到 summary.txt（若找不到返回 None）"""
    candidates = [
        os.path.join(summary_root, feature, "summary.txt"),
        os.path.join(summary_root, feature, "summary", "summary.txt"),
        os.path.join(summary_root, f"{feature}_summary.txt"),
        os.path.join(summary_root, f"{feature}.summary.txt"),
        os.path.join(summary_root, f"{feature}.txt"),
        os.path.join(summary_root, f"{feature}", f"{feature}_summary.txt"),
        os.path.join(summary_root, f"{feature}", "summary.txt"),
    ]
    for p in candidates:
        if os.path.isfile(p):
            return p
    return None

def load_feature_bin_priorities(feature, summary_root, default_n_bins=6):
    """
    返回 (priority_map, n_bins)
    priority_map: dict mapping original bin_number (int, 1-based) -> priority_rank (int, 1..n_bins, 高准确性 -> n_bins)
    """
    summary_file = find_summary_file_for_feature(feature, summary_root)
    if summary_file is None:
        # 没找到 summary，返回 identity 映射（1..default_n_bins）
        print(f"[WARN] summary for feature '{feature}' not found under {summary_root}. Using identity priority (bin->bin).")
        return ({i: i for i in range(1, default_n_bins+1)}, default_n_bins)

    # 读取 summary，宽松地按空白分割
    try:
        df = pd.read_csv(summary_file, sep=r'\s+', engine='python', dtype=str)
    except Exception as e:
        print(f"[WARN] Failed to read summary file {summary_file}: {e}. Using identity mapping.")
        return ({i: i for i in range(1, default_n_bins+1)}, default_n_bins)

    if df.shape[0] == 0:
        print(f"[WARN] Summary file {summary_file} is empty. Using identity mapping.")
        return ({i: i for i in range(1, default_n_bins+1)}, default_n_bins)

    # 找到 "Average" 行（第一列通常为 INV_ID）
    first_col = df.columns[0]
    mask_avg = df[first_col].astype(str).str.strip().str.lower() == 'average'
    if mask_avg.any():
        avg_row = df[mask_avg].iloc[0]
    else:
        # 回退：使用最后一行
        avg_row = df.iloc[-1]
    # 找到 bin 列（形如 bin_1 或 bin1）
    bin_cols = [c for c in df.columns if re.search(r'bin[_]?\d+$', c, flags=re.I)]
    if not bin_cols:
        # 如果没找到以 bin 开头的列，尝试匹配 '1','2',... 这种（不太常见）
        bin_cols = [c for c in df.columns if re.match(r'^\d+$', c)]
    if not bin_cols:
        print(f"[WARN] No bin columns found in {summary_file}. Using identity mapping.")
        return ({i: i for i in range(1, default_n_bins+1)}, default_n_bins)

    # 构建 bin -> accuracy 值
    acc_map = {}
    for c in bin_cols:
        m = re.search(r'(\d+)$', c)
        if not m:
            continue
        bin_idx = int(m.group(1))
        try:
            acc_val = float(avg_row[c])
        except Exception:
            # 如果解析失败，设为 0
            acc_val = 0.0
        acc_map[bin_idx] = acc_val

    # 确定 n_bins（以 acc_map 中最大 bin index 为准）
    n_bins = max(acc_map.keys()) if acc_map else default_n_bins

    # 对 bin 按 accuracy 降序排序，最高准确性的 bin 分配最高优先级 n_bins，次高 n_bins-1...
    sorted_bins = sorted(acc_map.items(), key=lambda x: (x[1], -x[0]), reverse=True)
    priority_map = {}
    # 如果少于 n_bins 的 bin 出现在 acc_map，仍需保证所有 1..n_bins 都有映射（未出现的设为 0）
    rank_val = n_bins
    for bin_idx, acc in sorted_bins:
        priority_map[bin_idx] = rank_val
        rank_val -= 1
    # 对于没有出现在 acc_map 的 bin，按原序补上（优先级从剩余较低值开始）
    for b in range(1, n_bins+1):
        if b not in priority_map:
            priority_map[b] = max(1, rank_val)
            rank_val -= 1

    return priority_map, n_bins

# ------------------------- 生成 SNP_BIN -------------------------
def process_position_bins(df, num_bins=1000):
    df['POS'] = df['SNP'].str.split('_').str[-1].astype(int)
    pos_values = df['POS'].values
    min_pos = pos_values.min()
    max_pos = pos_values.max()
    bins = np.linspace(min_pos, max_pos, num_bins + 1)

    # Add SNP_BIN column to the DataFrame
    df['SNP_BIN'] = pd.cut(df['POS'], bins=bins, labels=False, include_lowest=True) + 1

    mid_snp_map = {}
    for bin_id in range(1, num_bins + 1): # Iterate through 1-based bin_ids
        bin_start = bins[bin_id-1]
        bin_end = bins[bin_id]
        mask = (df['SNP_BIN'] == bin_id)
        sub_df = df[mask].copy()
        if not sub_df.empty:
            midpoint = (bin_start + bin_end) / 2
            sub_df['dist'] = (sub_df['POS'] - midpoint).abs()
            mid_snp = sub_df.nsmallest(1, 'dist').iloc[0]['SNP']
            mid_snp_map[bin_id] = mid_snp
    return df, bins, mid_snp_map # Return the modified df

# ------------------------- SNP选择 -------------------------
def select_snps_globally_weighted(df, mid_snp_map, weight_dict, summary_root, default_feature_bins=6, target_snp_count=1000):
    """
    Globally selects a target number of SNPs across the entire chromosome based on weighted features.
    """
    feature_meta = {}
    for feat, w in weight_dict.items():
        priority_map, n_bins = load_feature_bin_priorities(feat, summary_root, default_n_bins=default_feature_bins)
        feature_meta[feat] = {'weight': float(w), 'priority_map': priority_map, 'n_bins': n_bins}
        feat_bin_col = f"{feat}_BIN"
        if feat_bin_col not in df.columns:
            if feat in df.columns:
                target_bins = n_bins if feat.lower() != 'conserved' else 2
                try:
                    series_data = df[feat].fillna(0)
                    if feat.lower() == 'conserved':
                        df[feat_bin_col] = pd.cut(series_data, bins=target_bins, labels=False) + 1
                    else:
                        try:
                            qcut_res = pd.qcut(series_data, q=target_bins, labels=False, duplicates='drop')
                            if qcut_res.nunique() < target_bins:
                                df[feat_bin_col] = pd.cut(series_data, bins=target_bins, labels=False) + 1
                            else:
                                df[feat_bin_col] = qcut_res + 1
                        except ValueError:
                            df[feat_bin_col] = pd.cut(series_data, bins=target_bins, labels=False) + 1
                except Exception:
                    df[feat_bin_col] = 0
            else:
                df[feat_bin_col] = 0

    # Initialize weighted_score for all SNPs
    df['weighted_score'] = np.zeros(len(df))

    # Calculate weighted score for each SNP
    for feat, meta in feature_meta.items():
        w = meta['weight']
        pmap = meta['priority_map']
        feat_bin_col = f"{feat}_BIN"
        if feat_bin_col in df.columns:
            # Use .loc for safe assignment
            for bin_k, priority in pmap.items():
                if priority == 0:
                    continue
                mask = (df[feat_bin_col].fillna(0).astype(int) == int(bin_k))
                df.loc[mask, 'weighted_score'] += priority * w

    # Calculate distance to SNP_BIN midpoint for tie-breaking
    df['dist'] = np.nan
    for snp_bin, mid_snp_str in mid_snp_map.items():
        mid_pos = int(mid_snp_str.split('_')[-1])
        mask = (df['SNP_BIN'] == snp_bin)
        df.loc[mask, 'dist'] = (df.loc[mask, 'POS'] - mid_pos).abs()

    # Sort by weighted_score (descending) and then by dist (ascending)
    # Fill NaN distances with a large number so they are sorted last
    df['dist'] = df['dist'].fillna(df['dist'].max() + 1) # Ensure NaNs are at the end if any

    # Select top N SNPs
    selected_df = df.sort_values(by=['weighted_score', 'dist'], ascending=[False, True])
    selected_snps = selected_df['SNP'].head(target_snp_count).tolist()

    return selected_snps

# ------------------------- 文件输出 -------------------------
def save_snp_list(output_dir, snp_list):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    snp_file = output_dir / "selected_snps.txt"
    pd.Series(snp_list).to_csv(snp_file, index=False, header=False)
    print(f"Saved selected SNPs: {snp_file}")
    return snp_file

# ------------------------- Imputation Sub-routines -------------------------

def _run_plink_standardization(args, exp_dir, chrom, logger):
    # Create a directory for cached standardized benchmark files
    std_benchmark_dir = Path(args.output_dir) / "standardized_benchmark_cache"
    std_benchmark_dir.mkdir(exist_ok=True)
    std_benchmark_vcf = std_benchmark_dir / f"benchmark.chr{chrom}.standardized.vcf.gz"
    lock_path = std_benchmark_dir / f"benchmark.chr{chrom}.lock"

    # Use a file lock to prevent race conditions when creating the standardized benchmark
    lock = FileLock(lock_path)
    with lock:
        # Step 1: Create a standardized version of the benchmark VCF for the entire chromosome if it doesn't exist.
        # This is a one-time operation per chromosome.
        temp_prefix = std_benchmark_dir / f"temp_chr{chrom}"
        temp_vcf_path = temp_prefix.with_suffix('.vcf')
        try:
            if not std_benchmark_vcf.exists():
                logger.info(f"Creating standardized benchmark for chr{chrom} at {std_benchmark_vcf}...")
                pgen_prefix = temp_prefix.with_suffix('.pgen')
                plink_cmd_1_pgen = [
                    str(args.plink),
                    "--vcf", str(args.benchmark),
                    "--chr", str(chrom),
                    "--set-all-var-ids", "chr@_#",
                    "--new-id-max-allele-len", "23",
                    "--make-pgen",
                    "--out", str(temp_prefix),
                    "--max-alleles", "2",
                    "--autosome-num", "39",
                    "--allow-extra-chr"
                ]
                run_cmd(plink_cmd_1_pgen, logger)

                plink_cmd_1_export = [
                    str(args.plink),
                    "--pfile", str(temp_prefix),
                    "--export", "vcf", "bgz", # Add "bgz" to export directly to compressed VCF
                    "--out", str(temp_prefix),
                    "--autosome-num", "39",
                    "--allow-extra-chr"
                ]
                run_cmd(plink_cmd_1_export, logger)

                logger.info("Waiting for 5 seconds for the filesystem to sync...")
                time.sleep(5)
                # PLINK2 should now directly output .vcf.gz, so we just rename and index
                # Move the plink-generated .vcf.gz to std_benchmark_vcf
                plink_output_vcf_gz = temp_prefix.with_suffix('.vcf.gz')
                if plink_output_vcf_gz.exists():
                    os.rename(plink_output_vcf_gz, std_benchmark_vcf)
                else:
                    logger.error(f"PLINK2 did not generate expected output: {plink_output_vcf_gz}")
                    raise FileNotFoundError(f"PLINK2 output {plink_output_vcf_gz} not found.")

                run_cmd([str(args.bcftools), "index", "--threads", str(args.threads), str(std_benchmark_vcf)], logger, is_bcftools=True, bcftools_path=args.bcftools)
                logger.info(f"Standardized benchmark for chr{chrom} created successfully.")
            else:
                logger.info(f"Standardized benchmark for chr{chrom} already exists. Skipping generation.")
        finally:
            # Clean up temporary files from the first PLINK run
            for f in temp_prefix.parent.glob(f"{temp_prefix.name}.*"):
                if f.suffix not in ['.gz', '.csi'] and f.exists():
                    f.unlink()

    # Step 2: Use the chromosome-standardized benchmark VCF to extract SNPs for the current chunk.
    output_prefix = exp_dir / "standardized_chunk"
    standardized_vcf_gz = output_prefix.with_suffix(".vcf.gz")
    if not standardized_vcf_gz.exists():
        try:
            logger.info(f"Extracting SNPs for the current chunk from {std_benchmark_vcf}...")
            snp_file = exp_dir / "selected_snps.txt"
            
            # Generate position list for bcftools -R
            # snp_file contains lines like "chr1_12345"
            pos_file = output_prefix.with_suffix(".positions.txt")
            with open(snp_file, 'r') as f, open(pos_file, 'w') as pos_f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    # Expecting format like chr1_12345 or 1_12345
                    # Split by last underscore to separate chrom and pos
                    parts = line.rsplit('_', 1)
                    if len(parts) == 2:
                        chrom_part = parts[0]
                        pos_part = parts[1]
                        # Ensure chrom format matches VCF (usually "1", not "chr1" in this dataset)
                        # We strip 'chr' prefix to match the VCF which uses numbers (verified via tool check)
                        if chrom_part.lower().startswith('chr'):
                            chrom_part = chrom_part[3:]
                        
                        pos_f.write(f"{chrom_part}\t{pos_part}\n")

            # Use bcftools to extract directly from standardized benchmark
            run_cmd([str(args.bcftools), "view", "--threads", str(args.threads), "-R", str(pos_file), "-Oz", "-o", str(standardized_vcf_gz), str(std_benchmark_vcf)], logger, is_bcftools=True, bcftools_path=args.bcftools)
            run_cmd([str(args.bcftools), "index", "--threads", str(args.threads), str(standardized_vcf_gz)], logger, is_bcftools=True, bcftools_path=args.bcftools)
            logger.info("BCFtools standardization complete.")
        finally:
            # Clean up temporary position file
            if pos_file.exists():
                pos_file.unlink()
            # Clean up other temps if any (though we don't create plink temps anymore)
    else:
        logger.info(f"Found existing standardized chunk VCF: {standardized_vcf_gz}")
    return standardized_vcf_gz



def _run_phasing(args, exp_dir, target_vcf, logger):
    phased_vcf = exp_dir / "phased.vcf.gz"
    if not phased_vcf.exists():
        logger.info("Running phasing with Beagle...")
        run_cmd(["java", f"-Xmx{args.beagle_mem}g", "-jar", args.beagle_jar, f"gt={target_vcf}", f"out={exp_dir / 'phased'}", "impute=true", "seed=9823", f"nthreads={args.threads}"], logger)
        logger.info("Phasing complete.")
    else:
        logger.info(f"Found existing phased VCF: {phased_vcf}")
    return phased_vcf

def _run_conform_gt(args, exp_dir, phased_vcf, chrom, logger):
    conformed_vcf = exp_dir / "conformed.vcf.gz"
    if not conformed_vcf.exists():
        logger.info("Running Conform-GT...")
        run_cmd(["java", "-jar", args.conform_jar, f"ref={args.reference_vcf}", f"gt={phased_vcf}", f"chrom={chrom}", "match=POS", f"out={exp_dir / 'conformed'}"], logger)
    else:
        logger.info(f"Found existing conformed VCF: {conformed_vcf}")
    return conformed_vcf

def _run_beagle_imputation(args, exp_dir, conformed_vcf, chrom, start, end, logger):
    imputed_vcf = exp_dir / "imputed.vcf.gz"
    if not imputed_vcf.exists():
        logger.info("Running imputation with Beagle...")
        window_padding = 500000  # 500kb padding
        window_start = max(1, start - window_padding)
        window_end = end + window_padding
        beagle_cmd = [
            "java", f"-Xmx{args.beagle_mem}g", "-jar", args.beagle_jar,
            f"gt={conformed_vcf}",
            f"ref={args.reference_vcf}",
            f"chrom={chrom}:{window_start}-{window_end}",
            "impute=true",
            "ne=100",
            "seed=9823",
            f"out={exp_dir / 'imputed'}",
            f"nthreads={args.threads}",
        ]
        run_cmd(beagle_cmd, logger)
    else:
        logger.info(f"Found existing imputed VCF: {imputed_vcf}")
    return imputed_vcf

def _calculate_accuracy(args, exp_dir, imputed_vcf, original_benchmark_path, logger):
    # This function now expects to run on the merged VCF in the post-processing stage.
    # It generates results in exp_dir.
    
    # We construct the expected result file path based on the logic in imputation_accuracy.sh
    # script: OUTPUT_DIR/${out_name}_per_sample_results.txt
    imputed_name = imputed_vcf.name.replace('.vcf.gz', '')
    accuracy_result = exp_dir / f"{imputed_name}_per_sample_results.txt"
    
    def _is_valid_accuracy_file(p: Path) -> bool:
        try:
            if not p.exists():
                return False
            if p.stat().st_size <= 0:
                return False
            df = pd.read_csv(p, sep="\t")
            if df.shape[0] <= 0:
                return False
            required_cols = {"imputed_ids", "F-score"}
            if not required_cols.issubset(set(df.columns)):
                return False
            return True
        except Exception:
            return False

    if _is_valid_accuracy_file(accuracy_result):
        logger.info(f"[SKIP] Found existing accuracy results: {accuracy_result}")
        return True

    if accuracy_result.exists():
        try:
            logger.warning(f"Existing accuracy result file is invalid, will re-compute: {accuracy_result}")
            accuracy_result.unlink()
        except Exception as e:
            logger.warning(f"Failed to remove invalid accuracy result file {accuracy_result}: {e}")

    if not accuracy_result.exists():
        logger.info(f"Calculating accuracy for {imputed_vcf}...")
        
        # Determine the correct path for the accuracy script
        accuracy_script_path = Path(args.accuracy_script)
        if not accuracy_script_path.exists():
             # Fallback: check standard location if not found
             fallback_path = Path(os.environ.get("IDEA_ACCURACY_SCRIPT", "imputation_accuracy.sh"))
             if fallback_path.exists():
                 logger.warning(f"Accuracy script not found at {args.accuracy_script}, using fallback: {fallback_path}")
                 accuracy_script_path = fallback_path
             else:
                 logger.error(f"Accuracy script not found at {args.accuracy_script} or fallback location. Skipping accuracy calculation.")
                 return False

        # Ensure BCF exists for the benchmark
        # The script usually expects a BCF file or generates it. 
        # To avoid library errors in the shell script, we generate it here in Python if missing.
        bcf_gz_file = Path(str(original_benchmark_path).replace('.vcf.gz', '.bcf.gz'))
        
        if not bcf_gz_file.exists():
            logger.info(f"Generating BCF file from {original_benchmark_path}...")
            try:
                # Use bcftools to convert VCF to BCF
                # We use run_cmd so it handles LD_LIBRARY_PATH correctly
                run_cmd([
                    str(args.bcftools), "view", 
                    "-Ob", 
                    "-o", str(bcf_gz_file), 
                    "--threads", "20",
                    str(original_benchmark_path)
                ], logger, is_bcftools=True, bcftools_path=args.bcftools)
                
                # Index the BCF file
                run_cmd([
                    str(args.bcftools), "index", 
                    "--threads", "20", 
                    str(bcf_gz_file)
                ], logger, is_bcftools=True, bcftools_path=args.bcftools)
                logger.info(f"BCF file generated: {bcf_gz_file}")
            except Exception as e:
                logger.error(f"Failed to generate BCF file from {original_benchmark_path}: {e}")
                return False
        else:
            logger.info(f"Found existing BCF file: {bcf_gz_file}")

        # Run the accuracy script
        # We pass cwd=exp_dir so temporary files (like xac.vcf.gz) are created there, avoiding race conditions
        try:
            run_cmd([
                str(accuracy_script_path), 
                "-i", str(imputed_vcf), 
                "-w", str(original_benchmark_path),
                "-bw", str(bcf_gz_file), 
                "-t", "10", 
                "-o", f"{exp_dir}/"
            ], logger, is_bcftools=False, cwd=exp_dir, extra_env={
                "PYTHON_BIN": os.environ.get("IDEA_PYTHON_BIN", sys.executable),
                "BCFTOOLS_BIN": str(args.bcftools),
            })
        except Exception as e:
             logger.error(f"Accuracy script failed: {e}")
             return False
        
        if _is_valid_accuracy_file(accuracy_result):
            logger.info("[SUCCESS] Accuracy calculation complete.")
            return True

        logger.warning(f"[WARNING] Accuracy script finished but result file {accuracy_result} is missing or invalid.")
        return False

    return False

def _cleanup_intermediate_files(exp_dir, logger, final_outputs=None):
    logger.info("[INFO] Cleaning up intermediate files...")
    if final_outputs is None:
        final_outputs = []
    
    # Convert final_outputs to a set of absolute paths for efficient lookup
    final_outputs_set = {Path(f).resolve() for f in final_outputs}

    try:
        # List all files in the experiment directory
        all_files = list(exp_dir.iterdir())
        
        # Define patterns for files to always keep (e.g., logs, selected_snps.txt)
        always_keep_patterns = [
            "*.log",
            "selected_snps.txt",
            "*.positions.txt", # Keep position files for debugging if needed
            "*.vcf.gz.tbi", # Keep index files for final VCFs
            "*.vcf.gz.csi", # Keep index files for final VCFs
        ]
        
        for f_path in all_files:
            # Resolve path to handle symlinks and ensure comparison with absolute paths
            resolved_f_path = f_path.resolve()

            # Check if it's one of the explicitly final outputs
            if resolved_f_path in final_outputs_set:
                continue
            
            # Check if it matches any always_keep_patterns
            if any(f_path.match(pattern) for pattern in always_keep_patterns):
                continue

            # If it's not a final output and not an always_keep_pattern, and it's a file, delete it
            if f_path.is_file():
                try:
                    f_path.unlink()
                    logger.debug(f"Removed intermediate file: {f_path}")
                except Exception as e:
                    logger.warning(f"Failed to remove file {f_path}: {e}")
        logger.info("[SUCCESS] Intermediate files cleanup completed.")
    except Exception as e:
        logger.warning(f"[WARNING] Failed to clean some intermediate files: {e}")

# ------------------------- Imputation 流程 -------------------------
def run_imputation(args, exp_dir, snp_list, chrom, start, end, benchmark_vcf_path, logger, keep_intermediate=False):
    """Run the complete imputation pipeline for a set of selected SNPs."""
    exp_dir = Path(exp_dir)
    exp_dir.mkdir(parents=True, exist_ok=True)
    
    local_logger = logger.getChild(f'imputation_chr{chrom}_{exp_dir.name}')
    local_logger.info(f"Starting imputation in {exp_dir} for chromosome {chrom}")
    local_logger.debug(f"Beagle memory set to: {args.beagle_mem}g")
    local_logger.debug(f"Number of SNPs to process: {len(snp_list)}")

    try:
        snp_file = exp_dir / "selected_snps.txt"
        if not snp_file.exists():
            modified_snp_list = [f"chr{snp}" if not str(snp).startswith('chr') else str(snp) for snp in snp_list]
            pd.Series(modified_snp_list).to_csv(snp_file, index=False, header=False)
            local_logger.info(f"Saved SNPs to: {snp_file}")
        else:
            local_logger.info(f"File {snp_file} already exists, skipping SNP saving.")

        # Standardize the input VCF to get the genotypes for the selected SNPs for the current chunk
        target_vcf = _run_plink_standardization(args, exp_dir, chrom, local_logger)

        phased_vcf_path = exp_dir / "phased.vcf.gz"
        if phased_vcf_path.exists():
            local_logger.info(f"[SKIP] Phased VCF already exists: {phased_vcf_path}")
            phased_vcf = phased_vcf_path
        else:
            phased_vcf = _run_phasing(args, exp_dir, target_vcf, local_logger)

        final_vcf = phased_vcf # Initialize final_vcf with phased_vcf

        # Check if reference panel is provided for conform-gt and imputation
        if args.reference_vcf:
            conformed_vcf_path = exp_dir / "conformed.vcf.gz"
            if conformed_vcf_path.exists():
                local_logger.info(f"[SKIP] Conformed VCF already exists: {conformed_vcf_path}")
                conformed_vcf = conformed_vcf_path
            else:
                conformed_vcf = _run_conform_gt(args, exp_dir, phased_vcf, chrom, local_logger)
            final_vcf = conformed_vcf # Update final_vcf if conformed

            if args.skip_imputation:
                local_logger.info("Skipping imputation step as requested. Final VCF will be the conformed VCF.")
                # final_vcf is already conformed_vcf
            else:
                imputed_vcf_path = exp_dir / "imputed.vcf.gz"
                if imputed_vcf_path.exists():
                    local_logger.info(f"[SKIP] Imputed VCF already exists: {imputed_vcf_path}")
                    imputed_vcf = imputed_vcf_path
                else:
                    imputed_vcf = _run_beagle_imputation(args, exp_dir, conformed_vcf, chrom, start, end, local_logger)
                final_vcf = imputed_vcf
        else: # No reference panel provided
            if args.skip_imputation:
                local_logger.info("No reference panel provided and imputation is skipped. Final VCF will be the phased VCF.")
                # final_vcf is already phased_vcf
            else:
                local_logger.error("Reference panel (--reference-vcf) is required for imputation but not provided. Exiting.")
                raise ValueError("Reference panel missing for imputation.")

        if not final_vcf or not final_vcf.exists() or final_vcf.stat().st_size == 0:
            local_logger.error("Final VCF (imputed, conformed, or phased) is missing or empty. Processing failed.")
            return False

        # Index the final VCF file if it's not already indexed
        if not (Path(str(final_vcf) + ".tbi").exists() or Path(str(final_vcf) + ".csi").exists()):
            run_cmd([str(args.bcftools), "index", "--threads", "100", str(final_vcf)], local_logger, is_bcftools=True, bcftools_path=args.bcftools)

        # Ensure the final output is accessible as 'imputed.vcf.gz' for the merge step,
        # even if we skipped imputation.
        standard_output_name = exp_dir / "imputed.vcf.gz"
        if final_vcf.name != "imputed.vcf.gz":
            if standard_output_name.exists():
                standard_output_name.unlink()
            # Use symlink relative path if possible, or absolute
            try:
                standard_output_name.symlink_to(final_vcf.name)
                local_logger.info(f"Created symlink: {standard_output_name} -> {final_vcf.name}")
            except OSError:
                # Fallback to hard link or copy if symlink fails
                import shutil
                shutil.copy(final_vcf, standard_output_name)
                # Also copy index if exists
                for ext in ['.tbi', '.csi']:
                    src_idx = final_vcf.with_suffix(final_vcf.suffix + ext)
                    if src_idx.exists():
                        shutil.copy(src_idx, standard_output_name.with_suffix(standard_output_name.suffix + ext))
                local_logger.info(f"Copied final VCF to {standard_output_name}")

            # Ensure the symlinked/copied file is indexed
            if not (Path(str(standard_output_name) + ".tbi").exists() or Path(str(standard_output_name) + ".csi").exists()):
                 run_cmd([str(args.bcftools), "index", "--threads", "100", str(standard_output_name)], local_logger, is_bcftools=True, bcftools_path=args.bcftools)

        # Accuracy calculation is deferred to the post-processing stage (on the merged VCF).
        # We do NOT calculate accuracy per-chunk anymore.

        if not keep_intermediate:
            _cleanup_intermediate_files(exp_dir, local_logger, final_outputs=[str(final_vcf)])
        
        return True

    except subprocess.CalledProcessError as e:
        local_logger.error(f"[FATAL] A command failed in {exp_dir}: {str(e)}", exc_info=True)
        if hasattr(e, 'stdout') and e.stdout:
            local_logger.error(f"[FATAL] STDOUT: {e.stdout}")
        if hasattr(e, 'stderr') and e.stderr:
            local_logger.error(f"[FATAL] STDERR: {e.stderr}")
        return False
    except Exception as e:
        local_logger.error(f"[FATAL] An unexpected error occurred in {exp_dir}: {e}", exc_info=True)
        return False

# ------------------------- 主流程 -------------------------
def get_snp_count(bcftools_path, vcf_path, chrom):
    """使用bcftools获取VCF文件中的SNP数量"""
    try:
        cmd = [bcftools_path, 'index', '-n', vcf_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.strip().split('\n'):
            parts = line.split('\t')
            if parts[0] == str(chrom):
                return int(parts[2])
        return 0
    except Exception as e:
        print(f"[WARN] Could not determine SNP count for chromosome {chrom}: {e}")
        return 0

def get_chrom_length(bcftools_path, vcf_path, chrom):
    """使用bcftools获取染色体长度"""
    try:
        cmd = [bcftools_path, 'query', '-f', '%POS', '-r', str(chrom), vcf_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        positions = [int(p) for p in result.stdout.strip().split('\n') if p]
        return max(positions) if positions else 0
    except Exception as e:
        print(f"[WARN] Could not determine length for chromosome {chrom}: {e}")
        return 0

def process_chunk(args, chrom, start, end, benchmark_vcf_path, feature_df, weight_dict):
    """Processes a single genomic chunk."""
    # If chunking is disabled, use a simpler directory name (chrX) instead of chrX_start_end
    if args.enable_chunk:
        task_id = f"chr{chrom}_{start}_{end}"
    else:
        task_id = f"chr{chrom}"

    log_dir = Path(args.output_dir) / 'logs'
    log_dir.mkdir(exist_ok=True)
    logger = setup_logger(task_id, log_dir / f"{task_id}.log")
    logger.info(f"Processing region: {task_id}")

    try:
        # Filter feature_df to include only SNPs within the current chunk's genomic range
        # and that are present in the global_snp_list
        df_chunk = feature_df[
            (feature_df['CHR'] == chrom) &
            (feature_df['POS'] >= start) &
            (feature_df['POS'] < end) &
            (feature_df['SNP'].isin(args.global_snp_list)) # Filter by global_snp_list
        ].copy()

        if df_chunk.empty:
            logger.warning("No features found for this chunk after global SNP filtering.")
            return True

        # The following lines for SNP_BIN and select_snps_weighted are no longer needed
        # as global selection has already occurred.
        # We just need to ensure the snp_list passed to run_imputation is from df_chunk
        selected_snps = df_chunk['SNP'].tolist()

        if not selected_snps:
            logger.warning("No SNPs selected for this chunk after global SNP filtering.")
            return True

        chunk_output_dir = Path(args.output_dir) / task_id
        chunk_output_dir.mkdir(parents=True, exist_ok=True)

        success = run_imputation(
            args=args,
            exp_dir=chunk_output_dir,
            snp_list=selected_snps,
            chrom=chrom,
            start=start,
            end=end,
            benchmark_vcf_path=benchmark_vcf_path,
            logger=logger,
            keep_intermediate=args.keep_intermediate
        )

        if success:
            logger.info("Chunk processed successfully.")
        else:
            logger.warning("Chunk processing failed.")
        return success
    except Exception as e:
        logger.error(f"Error processing chunk {task_id}: {e}", exc_info=True)
        return False

def main():
    args = parse_args()
    if getattr(args, "select_only", False):
        run_select_only(args)
        return
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    main_log_file = output_dir / 'imputation_pipeline.log'
    logger = setup_logger('imputation_pipeline', main_log_file)

    logger.info("Starting Imputation Pipeline")
    logger.info(f"Arguments: {vars(args)}")

    try:
        pipeline = ImputationPipeline(args, logger)
        pipeline.run()
        logger.info("Imputation Pipeline finished successfully.")
    except Exception as e:
        logger.error("An unhandled exception occurred in the main pipeline.", exc_info=True)
        sys.exit(1)

class ImputationPipeline:
    def __init__(self, args, logger):
        self.args = args
        self.logger = logger
        self._check_dependencies()
        self.original_benchmark_path = Path(args.benchmark)
        self.preprocessed_vcfs = {}
        self.feature_df = load_features(self.args.feature_file)
        self.weight_dict = load_feature_weights(self.args.feature_weight)
        self.selection_chromosomes = self._parse_selection_chromosomes()
        self.test_chromosome = self._parse_test_chromosome()
        self.chromosomes = [self.test_chromosome] if self.test_chromosome is not None else list(self.selection_chromosomes)
        self.selected_snps_by_chrom = {}
        self.selected_snps_genome = []
        self.running_tasks = {} # future -> (task, retries, mem_gb)
        self.current_mem_gb = 0

    def _check_dependencies(self):
        """Checks if all required external tools are available."""
        self.logger.info("--- Checking dependencies ---")
        dependencies = {
            "bcftools": self.args.bcftools,
            "plink": self.args.plink,
            "beagle.jar": self.args.beagle_jar
        }
        all_found = True
        for name, path in dependencies.items():
            if not Path(path).exists():
                self.logger.error(f"Dependency not found: {name} at path '{path}'")
                all_found = False
            else:
                self.logger.info(f"✅ Found {name} at {path}")

        if not all_found:
            self.logger.error("One or more dependencies are missing. Please check the paths and try again.")
            sys.exit(1)
        self.logger.info("--- All dependencies found ---")

    def run(self):
        self.logger.info("Pipeline run started.")
        self._select_genome_snps()
        self._preprocess_vcfs()
        all_tasks = self._generate_tasks()
        
        if not all_tasks:
            self.logger.info("No tasks to run.")
            results = {}
        else:
            results = self._process_tasks(all_tasks)

        self._post_process(results)
        self.logger.info("Pipeline run finished.")

    def _parse_selection_chromosomes(self):
        chromosomes = _parse_chr_arg(self.args.chr, feature_df=self.feature_df)
        if not chromosomes:
            self.logger.error("No chromosomes inferred from feature file. Please provide --chr.")
            sys.exit(1)
        self.logger.info(f"Chromosomes for SNP selection: {chromosomes}")
        return chromosomes

    def _parse_test_chromosome(self):
        test_chrom = _choose_test_chromosome(self.args.test_chr, self.feature_df, self.selection_chromosomes)
        if test_chrom is None:
            self.logger.info("Test chromosome: all selected chromosomes (no --test-chr).")
            return None
        if int(test_chrom) not in set(int(x) for x in self.selection_chromosomes):
            self.logger.error(f"--test-chr={test_chrom} is not included in selection chromosomes {self.selection_chromosomes}.")
            sys.exit(1)
        self.logger.info(f"Test chromosome: {int(test_chrom)}")
        return int(test_chrom)

    def _select_genome_snps(self):
        self.logger.info("--- SNP selection (genome design) ---")
        chrom_to_n = {}
        for chrom in self.selection_chromosomes:
            chrom_to_n[int(chrom)] = int((self.feature_df['CHR'] == chrom).sum())

        per_chrom_targets = {}
        if getattr(self.args, "target_snp_total", None) is not None:
            per_chrom_targets = _allocate_target_counts(chrom_to_n, self.args.target_snp_total)
            self.logger.info(f"Target SNP selection mode: total={int(self.args.target_snp_total)} across chromosomes -> per-chrom targets: {per_chrom_targets}")
        else:
            self.logger.info(f"Target SNP selection mode: per-chrom={int(self.args.target_snp_count)}")

        selected_all = []
        for chrom in self.selection_chromosomes:
            chrom_df = self.feature_df[self.feature_df['CHR'] == chrom].copy()
            if chrom_df.empty:
                continue
            target_for_chrom = int(per_chrom_targets.get(int(chrom), self.args.target_snp_count)) if per_chrom_targets else int(self.args.target_snp_count)
            if target_for_chrom <= 0:
                self.selected_snps_by_chrom[int(chrom)] = []
                continue
            snps = _select_snps_stratified_by_windows(
                chrom_df=chrom_df,
                weight_dict=self.weight_dict,
                summary_root=self.args.summary_root,
                default_feature_bins=self.args.default_feature_bins,
                window_bp=self.args.snp_bin_size,
                target_snp_count=target_for_chrom,
            )
            self.selected_snps_by_chrom[int(chrom)] = snps
            selected_all.extend(snps)

        self.selected_snps_genome = selected_all
        out_dir = Path(self.args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        genome_snp_file = out_dir / "selected_snps.genome.txt"
        pd.Series(self.selected_snps_genome).to_csv(genome_snp_file, index=False, header=False)
        self.logger.info(f"Saved genome SNP list: {genome_snp_file} (n={len(self.selected_snps_genome)})")

    def _preprocess_vcfs(self):
        preprocessed_dir = Path(self.args.output_dir) / "preprocessed_vcfs"
        preprocessed_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info("--- Pre-processing: Splitting VCFs by chromosome ---")

        with ThreadPoolExecutor(max_workers=self.args.threads) as executor:
            futures = {executor.submit(self._process_chromosome_vcf, chrom, preprocessed_dir): chrom for chrom in self.chromosomes}
            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    ref_vcf_out, bench_vcf_out = future.result()
                    self.preprocessed_vcfs[chrom] = {'ref': str(ref_vcf_out), 'bench': str(bench_vcf_out)}
                except Exception as exc:
                    self.logger.error(f'Chromosome {chrom} generated an exception during VCF preprocessing: {exc}')

        self.logger.info("--- VCF pre-processing complete ---")

    def _process_chromosome_vcf(self, chrom, preprocessed_dir):
        self.logger.info(f"Processing chromosome {chrom}...")
        ref_vcf_out = preprocessed_dir / f"ref_chr{chrom}.vcf.gz"

        # This logic assumes input paths are templates like '/path/to/file_chr{chrom}.vcf.gz'
        # and creates symlinks instead of extracting from a monolithic VCF.
        def link_file_and_index(vcf_type, source_template, dest_path):
            try:
                source_path = Path(str(source_template).format(chrom=chrom))
            except KeyError:
                self.logger.error(f"The {vcf_type} VCF path '{source_template}' does not seem to be a template with '{{chrom}}'. Cannot proceed with linking.")
                raise

            if not dest_path.exists():
                self.logger.info(f"Linking {vcf_type} VCF for chr {chrom}: {source_path} -> {dest_path}")
                if not source_path.exists():
                    raise FileNotFoundError(f"Source {vcf_type} VCF not found: {source_path}")
                
                # Create symlink for the main VCF file
                dest_path.symlink_to(source_path.resolve())

                # Attempt to link the index file (.tbi or .csi)
                found_index = False
                for ext in [".tbi", ".csi"]:
                    source_idx = source_path.with_suffix(source_path.suffix + ext)
                    if source_idx.exists():
                        dest_idx = dest_path.with_suffix(dest_path.suffix + ext)
                        dest_idx.symlink_to(source_idx.resolve())
                        self.logger.info(f"Linked index file {source_idx.name}")
                        found_index = True
                        break
                if not found_index:
                    self.logger.warning(f"Index file for {source_path} not found. Downstream tools may need to create one.")

        bench_vcf_out = preprocessed_dir / f"bench_chr{chrom}.vcf.gz"

        link_file_and_index("reference", self.args.reference_vcf, ref_vcf_out)
        link_file_and_index("benchmark", self.original_benchmark_path, bench_vcf_out)

        return ref_vcf_out, bench_vcf_out

    def _generate_tasks(self):
        tasks = []
        for chrom in self.chromosomes:
            ref_vcf = self.preprocessed_vcfs[chrom]['ref']
            
            # Unified chunking logic
            max_snps = self.args.max_snps_per_chunk
            max_size = self.args.max_chunk_size if self.args.dynamic_window else self.args.window_size
            mode = "dynamic (size/count)" if self.args.dynamic_window else "fixed (window size)"

            self.logger.info(f"Generating chunks for chr {chrom} in {mode} mode with max_snps={max_snps} and max_size={max_size}bp.")

            # Filter feature_df for the current chromosome
            chrom_feature_df = self.feature_df[self.feature_df['CHR'] == chrom].copy()

            if chrom_feature_df.empty:
                self.logger.warning(f"No features found for chromosome {chrom}. Skipping global SNP selection.")
                continue

            # The feature universe defines the selectable genomic range and is
            # already loaded in memory. Using its maximum coordinate avoids an
            # unnecessary scan of millions of reference records. For the IDEA
            # contract this is the correct boundary: no marker can be selected
            # outside the feature universe.
            chrom_start = 1
            chrom_end = int(chrom_feature_df['POS'].max())
            self.logger.info(
                f"Chromosome {chrom} range: {chrom_start}-{chrom_end} "
                "(from feature universe)."
            )

            global_snp_list = self.selected_snps_by_chrom.get(int(chrom), [])

            if not global_snp_list:
                self.logger.warning(f"No SNPs selected globally for chromosome {chrom}. Skipping chunk generation.")
                continue
            
            # Create a sorted list of selected SNP positions for quick lookup
            # This is used to ensure we have enough "feature" SNPs in a chunk if needed, 
            # though the primary driver now is contiguous coverage.
            # Actually, the original logic used global_snp_list to drive chunking.
            # Now we use physical/count limits but ensure contiguous blocks.
            
            # To efficiently find SNP counts, we can use bisect on the full 'positions' list (all variants in Ref)
            # OR 'global_snp_positions' (selected variants).
            # The 'max_snps' usually refers to the number of variants processed in imputation (Ref panel density),
            # or the number of target SNPs? Usually it's about computational load, so Ref variants count matters most.
            # However, previous logic used 'global_snp_positions' (selected features).
            # Let's stick to using 'positions' (all Ref variants) for density control if that's what matters for Beagle,
            # OR 'global_snp_positions' if we only care about feature density.
            # Given this is "imputation", Beagle performance depends on total markers. 
            # But let's look at the old code: it used `global_snp_positions`.
            # Let's stick to `global_snp_positions` for the "count" limit to minimize behavior change,
            # BUT we must use `chrom_start` to `chrom_end` for coverage.
            
            task_args = argparse.Namespace(**vars(self.args))
            task_args.reference_vcf = ref_vcf
            task_args.global_snp_list = global_snp_list
            benchmark_vcf_path = self.preprocessed_vcfs[chrom]['bench']

            if not self.args.enable_chunk:
                self.logger.info(f"Chunking disabled (default). Generating single task for chromosome {chrom} ({chrom_start}-{chrom_end}).")
                tasks.append((task_args, chrom, chrom_start, chrom_end, benchmark_vcf_path, self.feature_df, self.weight_dict))
            else:
                self.logger.info(f"Chunking enabled. Generating chunks for chromosome {chrom}.")
                global_snp_positions = sorted(chrom_feature_df[chrom_feature_df['SNP'].isin(global_snp_list)]['POS'].tolist())
                
                import bisect

                current_chunk_start = chrom_start
                
                while current_chunk_start <= chrom_end:
                    # 1. Determine end by size
                    candidate_end_size = current_chunk_start + max_size
                    
                    # 2. Determine end by SNP count (using selected SNPs as proxy for density/feature richness)
                    # Find the index of the first SNP >= current_chunk_start
                    start_snp_idx = bisect.bisect_left(global_snp_positions, current_chunk_start)
                    
                    # We want a chunk to contain at most max_snps
                    target_snp_idx = start_snp_idx + max_snps
                    
                    if target_snp_idx < len(global_snp_positions):
                        # The position of the (max_snps + 1)-th SNP
                        # We set the end BEFORE this SNP to keep count <= max_snps
                        # Actually, let's just use the position of that SNP as the limit.
                        candidate_end_count = global_snp_positions[target_snp_idx]
                    else:
                        candidate_end_count = chrom_end + 1 # Beyond end

                    # 3. Choose the smaller end, but don't go beyond chrom_end
                    chunk_end_pos = min(candidate_end_size, candidate_end_count)
                    
                    # Cap at chrom_end
                    if chunk_end_pos > chrom_end:
                        chunk_end_pos = chrom_end

                    # Ensure we make progress (at least 1 bp, or jump to next SNP if empty?)
                    if chunk_end_pos <= current_chunk_start:
                         chunk_end_pos = current_chunk_start + 1 # Force at least 1bp progress

                    # 4. Create Task
                    # task_args and benchmark_vcf_path are already defined above
                    
                    tasks.append((task_args, chrom, current_chunk_start, chunk_end_pos, benchmark_vcf_path, self.feature_df, self.weight_dict))
                    
                    # 5. Prepare for next chunk (contiguous)
                    current_chunk_start = chunk_end_pos + 1

        self.logger.info(f"Generated {len(tasks)} tasks to be processed.")
        return tasks



    def _process_tasks(self, tasks):
        self.logger.info(f"Starting parallel imputation for {len(tasks)} genomic windows.")
        self.logger.info(f"Resource limits: Max concurrent tasks = {self.args.max_concurrent_tasks}, Max total memory = {self.args.max_total_mem} GB.")
        
        results = {}
        task_queue = list(tasks)
        
        with tqdm(total=len(tasks), desc="Processing Chunks") as pbar, \
             ThreadPoolExecutor(max_workers=self.args.threads) as executor:
            while task_queue or self.running_tasks:
                # Submit new tasks if resources are available
                while task_queue and \
                      len(self.running_tasks) < self.args.max_concurrent_tasks and \
                      self.current_mem_gb + self.args.beagle_mem <= self.args.max_total_mem:
                    
                    task_item = task_queue.pop(0)

                    # Unpack if it's a retried task: (original_task_tuple, retry_count)
                    if isinstance(task_item, tuple) and len(task_item) == 2 and isinstance(task_item[0], tuple):
                        task, retries_so_far = task_item
                    else:
                        # It's a new task
                        task = task_item
                        retries_so_far = 0

                    # Now 'task' should be the 7-element tuple of arguments for process_chunk
                    if not (isinstance(task, tuple) and len(task) == 7):
                        self.logger.error(f"Skipping malformed task item from queue: {task_item}")
                        continue

                    mem_for_task = self.args.beagle_mem

                    self.logger.info(f"Submitting task {task[1:4]} with {mem_for_task} GB memory. "
                                     f"Current tasks: {len(self.running_tasks)}, Current memory: {self.current_mem_gb} GB.")

                    future = executor.submit(process_chunk, *task)
                    # Store the task with its current retry count
                    self.running_tasks[future] = (task, retries_so_far, mem_for_task)
                    self.current_mem_gb += mem_for_task

                if not self.running_tasks:
                    time.sleep(1)
                    continue

                done, _ = wait(self.running_tasks.keys(), return_when=FIRST_COMPLETED)

                for future in done:
                    task, retries, mem_gb = self.running_tasks.pop(future)
                    self.current_mem_gb -= mem_gb
                    task_id = (task[1], task[2], task[3])

                    try:
                        success = future.result()
                        if success:
                            self.logger.info(f"✅ Task {task_id} completed successfully.")
                            results[task_id] = True
                            pbar.update(1)
                        else:
                            self.logger.warning(f"⚠️ Task {task_id} failed.")
                            if retries < self.args.max_retries:
                                self.logger.info(f"Retrying task {task_id} ({retries + 1}/{self.args.max_retries}). Adding back to queue.")
                                task_queue.insert(0, (task, retries + 1))
                            else:
                                self.logger.error(f"❌ Task {task_id} failed after {self.args.max_retries} retries.")
                                results[task_id] = False
                    except Exception as exc:
                        self.logger.error(f"❌ Task {task_id} generated an exception: {exc}", exc_info=True)
                        if retries < self.args.max_retries:
                            self.logger.info(f"Retrying task {task_id} ({retries + 1}/{self.args.max_retries}) due to exception. Adding back to queue.")
                            task_queue.insert(0, (task, retries + 1))
                        else:
                            self.logger.error(f"❌ Task {task_id} failed after {self.args.max_retries} retries due to exception.")
                            results[task_id] = False
                


        self.logger.info("All tasks have been processed.")
        return results

    def _post_process(self, results):
        self.logger.info("--- Post-processing results ---")
        
        # First, merge VCF chunks for each chromosome (ONLY if chunking was enabled)
        # If chunking was disabled, we still need to ensure the final VCF is in the right place/named correctly
        # for the accuracy check.
        
        if self.args.enable_chunk:
            merge_vcf_chunks(self.args.output_dir, self.args.bcftools, self.chromosomes, self.logger)
            if self.args.cleanup:
                cleanup_chunk_directories(self.args.output_dir, self.chromosomes, self.logger)
        else:
             # If chunking disabled, we need to locate the single output file and rename/link it 
             # to merged_imputed_chrX.vcf.gz so that accuracy check finds it.
             for chrom in self.chromosomes:
                 # Find the task output directory
                 # If chunking was disabled, the directory is just chr{chrom}
                 
                 single_task_dir = Path(self.args.output_dir) / f"chr{chrom}"
                 if not single_task_dir.exists():
                     # Fallback to old pattern if not found (though new runs should use chr{chrom})
                     chunk_dirs = glob.glob(str(Path(self.args.output_dir) / f"chr{chrom}_*_*"))
                     if chunk_dirs:
                         single_task_dir = Path(chunk_dirs[0])
                     else:
                         self.logger.warning(f"No output directory found for chromosome {chrom}.")
                         continue
                 
                 imputed_vcf = single_task_dir / "imputed.vcf.gz"
                 
                 if imputed_vcf.exists():
                     target_merged_vcf = Path(self.args.output_dir) / f"merged_imputed_chr{chrom}.vcf.gz"
                     try:
                         # Symlink or copy
                         if target_merged_vcf.exists() or target_merged_vcf.is_symlink():
                             target_merged_vcf.unlink()
                         target_merged_vcf.symlink_to(imputed_vcf.resolve())
                         self.logger.info(f"Linked single chunk result to {target_merged_vcf}")
                         
                         # Also link index
                         for ext in ['.tbi', '.csi']:
                             src_idx = imputed_vcf.with_suffix(imputed_vcf.suffix + ext)
                             if src_idx.exists():
                                 target_idx = target_merged_vcf.with_suffix(target_merged_vcf.suffix + ext)
                                 if target_idx.exists() or target_idx.is_symlink():
                                     target_idx.unlink()
                                 target_idx.symlink_to(src_idx.resolve())
                     except Exception as e:
                         self.logger.error(f"Failed to link result for chr {chrom}: {e}")

        # Then, run accuracy calculation on the merged VCFs
        self._run_final_accuracy_check()

        successful_tasks = all(results.values())
        if successful_tasks:
            self.logger.info("🎉 All windows processed successfully!")
        else:
            failed_windows = [task_id for task_id, success in results.items() if not success]
            self.logger.warning(f"⚠️ Some windows failed processing: {failed_windows}")

    def _run_final_accuracy_check(self):
        self.logger.info("--- Running Final Accuracy Calculation ---")
        
        # Check if accuracy script exists
        accuracy_script_path = Path(self.args.accuracy_script)
        if not accuracy_script_path.exists():
             # Fallback: check standard location if not found
             fallback_path = Path(os.environ.get("IDEA_ACCURACY_SCRIPT", "imputation_accuracy.sh"))
             if fallback_path.exists():
                 self.logger.warning(f"Accuracy script not found at {self.args.accuracy_script}, using fallback: {fallback_path}")
                 accuracy_script_path = fallback_path
             else:
                 self.logger.warning(f"Accuracy script not found at {self.args.accuracy_script}. Skipping accuracy calculation.")
                 return

        accuracy_dir = Path(self.args.output_dir) / "accuracy_results"
        accuracy_dir.mkdir(exist_ok=True)

        all_accuracy_dfs = []

        for chrom in self.chromosomes:
            merged_vcf = Path(self.args.output_dir) / f"merged_imputed_chr{chrom}.vcf.gz"
            if not merged_vcf.exists():
                self.logger.warning(f"Merged VCF for chromosome {chrom} not found. Skipping accuracy calculation.")
                continue

            # The original benchmark path for this chromosome
            original_benchmark_path = self.preprocessed_vcfs[chrom]['bench']
            
            # Use the directory containing the merged VCF as the output directory for accuracy results
            # merged_vcf is in self.args.output_dir (root output dir)
            chrom_acc_dir = Path(self.args.output_dir)

            self.logger.info(f"Calculating accuracy for merged VCF: {merged_vcf}")
            
            # Call the accuracy calculation helper
            # IMPORTANT: Pass the corrected accuracy_script_path string to the helper function
            # We temporarily update args.accuracy_script so the helper uses the correct path
            original_script_path = self.args.accuracy_script
            self.args.accuracy_script = str(accuracy_script_path)
            
            _calculate_accuracy(self.args, chrom_acc_dir, merged_vcf, original_benchmark_path, self.logger)
            
            # Restore original path just in case
            self.args.accuracy_script = original_script_path
            
            # Collect results for final summary report
            # The script output name depends on the script itself, usually derived from input name
            # Our script typically outputs: {OUTPUT_DIR}/{IMPUTED_PREFIX}_per_sample_results.txt
            
            imputed_prefix = merged_vcf.name.replace('.vcf.gz', '')
            result_file = chrom_acc_dir / f"{imputed_prefix}_per_sample_results.txt"
            
            if result_file.exists():
                try:
                    df = pd.read_csv(result_file, sep="\t")
                    df['chromosome'] = chrom
                    all_accuracy_dfs.append(df)
                    self.logger.info(f"Collected accuracy results for chr {chrom}")
                except Exception as e:
                    self.logger.warning(f"Could not read or process accuracy results for chr {chrom}: {e}")
            else:
                 self.logger.warning(f"Expected accuracy result file not found: {result_file}")

        if all_accuracy_dfs:
            final_df = pd.concat(all_accuracy_dfs, ignore_index=True)
            final_report_path = Path(self.args.output_dir) / "final_imputation_accuracy_report.csv"
            final_df.to_csv(final_report_path, index=False)
            self.logger.info(f"✅ Final accuracy report saved to: {final_report_path}")

def merge_vcf_chunks(output_dir, bcftools_path, chromosomes, logger):
    logger.info("--- Merging VCF chunks ---")
    for chrom in chromosomes:
        # Find all imputed.vcf.gz files in the subdirectories
        chunk_dirs_glob = str(Path(output_dir) / f"chr{chrom}_*_*" / "imputed.vcf.gz")
        chunk_vcfs = glob.glob(chunk_dirs_glob)

        if not chunk_vcfs:
            logger.info(f"No VCF chunks found for chromosome {chrom} using pattern: {chunk_dirs_glob}. Skipping.")
            continue

        # Sort the VCF files numerically/version-wise to ensure correct order
        chunk_vcfs.sort(key=lambda f: [int(x) for x in re.findall(r'\d+', Path(f).parent.name)])

        merged_vcf_path = Path(output_dir) / f"merged_imputed_chr{chrom}.vcf.gz"
        vcf_list_file = Path(output_dir) / f"chr{chrom}_vcf_list.txt"

        with open(vcf_list_file, 'w') as f:
            for vcf in chunk_vcfs:
                f.write(f"{vcf}\n")
        logger.info(f"Created VCF list file for merging: {vcf_list_file}")

        # Use bcftools concat for stitching chunks together (same samples, different sites)
        # --naive is faster if headers are compatible (usually true for chunks from same run)
        # But standard concat is safer. -a allows overlaps if any (though we shouldn't have overlaps ideally)
        # -D removes duplicates if any.
        
        # NOTE: We previously used 'merge' which is wrong for chunks of same samples. 'concat' is correct.
        merge_cmd = [
            str(bcftools_path), "concat",
            "--file-list", str(vcf_list_file),
            "-Oz", "-o", str(merged_vcf_path), "--threads", "100", "-a", "-D"
        ]
        logger.info(f"Running concat command for chromosome {chrom}: {' '.join(merge_cmd)}")
        try:
            run_cmd(merge_cmd, logger, is_bcftools=True, bcftools_path=bcftools_path)
            logger.info(f"Successfully merged VCFs into {merged_vcf_path}")

            # Create index for the merged file
            index_cmd = [str(bcftools_path), "index", "-t", str(merged_vcf_path), "--threads", "100"]
            logger.info(f"Indexing merged VCF for chromosome {chrom}: {' '.join(index_cmd)}")
            run_cmd(index_cmd, logger, is_bcftools=True, bcftools_path=bcftools_path)
            logger.info(f"Indexed {merged_vcf_path}")

        except subprocess.CalledProcessError as e:
            logger.error(f"Error merging VCFs for chromosome {chrom}: {e}")
            # Continue to the next chromosome
        finally:
            # Clean up VCF list file
            if os.path.exists(vcf_list_file):
                os.remove(vcf_list_file)
                logger.info(f"Removed temporary VCF list file: {vcf_list_file}")

    logger.info("--- VCF merging completed for all chromosomes ---")



def cleanup_chunk_directories(output_dir, chromosomes, logger):
    logger.info("--- Cleaning up chunk directories ---")
    import shutil
    for chrom in chromosomes:
        chunk_dirs = glob.glob(str(Path(output_dir) / f"chr{chrom}_*_*_*"))
        for chunk_dir in chunk_dirs:
            try:
                # 递归删除目录及其内容
                shutil.rmtree(chunk_dir)
                logger.info(f"Removed directory: {chunk_dir}")
            except Exception as e:
                logger.warning(f"Failed to remove directory {chunk_dir}: {e}")
    logger.info("✅ Cleanup of chunk directories completed.")

if __name__ == "__main__":
    main()
