#!/usr/bin/env python3
"""
generate_feature.py

生成 minimal feature 文件，列为：
CHR  SNP  A1  A2  MAF  Conserved  L2  Conserved_score

用法示例：
python3 generate_feature.py \
  --bfile examples/plink/TP_chr6 \
  --out results/TP_feature_minimal.txt \
  --conserved-bed /path/to/sus_scrofa_conserved.bed \
  --conscore-bed /path/to/sus_scrofa_conserved_score.bed \
  --plink plink \
  --ldsc /path/to/ldsc.py \
  --bedtools bedtools

说明：
- 如果已存在 `${prefix}_maf.frq`、`${prefix}_l2.txt` 或 `${prefix}_ld.l2.ldscore.gz`、以及 `${breed}_chr6_conscore.txt`，脚本会直接使用它们。
- 否则需提供对应工具路径 (--plink, --ldsc, --bedtools) 来生成中间文件。
"""

import argparse
import os
import sys
import gzip
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np

def run(cmd):
    print("[CMD] " + " ".join(cmd) if isinstance(cmd, list) else "[CMD] " + cmd)
    subprocess.run(cmd, check=True)

def ensure_maf(prefix, outdir, plink_bin):
    # expected out: {outdir}/{sample}_maf.frq  where sample is basename(prefix)
    sample = Path(prefix).name
    maf_path = Path(outdir) / f"{sample}_maf.frq"
    if maf_path.exists():
        print(f"[INFO] Found existing MAF: {maf_path}")
        return maf_path
    if plink_bin is None:
        print(f"[WARN] MAF file {maf_path} not found and --plink not provided. MAF will be filled with 0.")
        return None
    # run plink --bfile prefix --freq --out outdir/sample_maf
    out_prefix = str(Path(outdir) / f"{sample}_maf")
    cmd = [plink_bin, "--bfile", prefix, "--freq", "--out", out_prefix]
    run(cmd)
    if maf_path.exists():
        print(f"[INFO] Generated MAF: {maf_path}")
        return maf_path
    else:
        print(f"[WARN] Expected maf not produced: {maf_path}")
        return None

def read_maf(maf_path, bim_df):
    # plink .frq typically whitespace-delimited with header: CHR SNP A1 A2 MAF NCHROBS
    if maf_path is None:
        # return zeros
        return pd.Series(0.0, index=bim_df['SNP'])
    try:
        maf_df = pd.read_csv(maf_path, delim_whitespace=True, dtype={'SNP':str})
        maf_map = maf_df.set_index('SNP')['MAF'].to_dict()
    except Exception as e:
        print(f"[WARN] Failed to read maf file {maf_path}: {e}")
        maf_map = {}
    return pd.Series([maf_map.get(snp, 0.0) for snp in bim_df['SNP']], index=bim_df['SNP'])

def ensure_ldsc(prefix, outdir, ldsc_py):
    sample = Path(prefix).name
    ld_gz = Path(outdir) / f"{sample}_ld.l2.ldscore.gz"
    l2_txt = Path(outdir) / f"{sample}_l2.txt"
    if l2_txt.exists():
        print(f"[INFO] Found existing L2 mapping: {l2_txt}")
        return l2_txt
    if ld_gz.exists():
        print(f"[INFO] Found existing ldscore gz: {ld_gz}")
        # will create mapping below
        pass
    elif ldsc_py is None:
        print(f"[WARN] ldscore not found and --ldsc not provided. L2 will be filled with 0.")
        return None
    else:
        # run ldsc: python ldsc.py --bfile prefix --l2 --ld-wind-cm 1 --maf 0 --out outdir/sample_ld --yes-really
        out_pref = str(Path(outdir) / f"{sample}_ld")
        cmd = ["python", ldsc_py, "--bfile", prefix, "--l2", "--ld-wind-cm", "1", "--maf", "0", "--out", out_pref, "--yes-really"]
        run(cmd)
    # now ld_gz should exist
    if not ld_gz.exists():
        print(f"[WARN] ldscore .gz still not found: {ld_gz}")
        return None
    # parse ld_gz: second column variant id, last column ldscore
    ldmap = {}
    with gzip.open(ld_gz, 'rt') as fh:
        header = fh.readline()
        for ln in fh:
            parts = ln.strip().split()
            if len(parts) >= 2:
                var = parts[1]
                try:
                    val = float(parts[-1])
                except:
                    val = 0.0
                ldmap[var] = val
    # We will map to bim outside, return a dict
    return ldmap

def ensure_l2_txt(prefix, outdir, ldsc_py):
    # wrapper: prefer existing l2.txt file if present, else produce ldscore.gz and return map/dict
    sample = Path(prefix).name
    l2_txt = Path(outdir) / f"{sample}_l2.txt"
    if l2_txt.exists():
        print(f"[INFO] Found existing L2 txt: {l2_txt}")
        # read file into dict
        df = pd.read_csv(l2_txt, sep="\t", header=None, names=['SNP','L2'])
        return dict(zip(df['SNP'].astype(str), df['L2'].astype(float)))
    # else try produce/read ld.l2.ldscore.gz
    ldmap = ensure_ldsc(prefix, outdir, ldsc_py)
    if isinstance(ldmap, dict):
        return ldmap
    return None

def ensure_conserved_scores(prefix, outdir, breed, conserved_bed, conscore_bed, bedtools_bin):
    # create tss bed from bim (POS-1, POS) -> tss_file
    sample = Path(prefix).name
    bim = Path(prefix + ".bim")
    if not bim.exists():
        raise FileNotFoundError(f"BIM file not found: {bim}")
    tss_file = Path(outdir) / f"{breed}_chr6_tss.bed"
    if not tss_file.exists():
        # read bim: chrom, snp, cm, pos, a1, a2
        bim_df = pd.read_csv(bim, sep='\t', header=None, usecols=[0,1,3,4,5], names=['CHR','SNP','POS','A1','A2'])
        # POS is column 3? above uses 0-indexed usecols, careful: we used cols [0,1,3,4,5], so POS at col index 3
        # create START = POS-1, END = POS
        if 'POS' not in bim_df.columns:
            # try reading column 3 explicitly
            bim_df_full = pd.read_csv(bim, sep='\t', header=None)
            # assume POS at col index 3
            bim_df = bim_df_full.loc[:, [0,1,3,4,5]]
            bim_df.columns = ['CHR','SNP','POS','A1','A2']
        bim_df['START'] = bim_df['POS'] - 1
        bim_df['END'] = bim_df['POS']
        tss_df = bim_df[['CHR','START','END','SNP']]
        tss_df.to_csv(tss_file, sep='\t', header=False, index=False)
        print(f"[INFO] Created TSS bed: {tss_file}")
    # conserved score: prefer provided conscore file if exists under outdir
    conscore_txt = Path(outdir) / f"{breed}_chr6_conscore.txt"
    if conscore_txt.exists():
        print(f"[INFO] Found existing conscore: {conscore_txt}")
        # read later
        return conscore_txt
    # if conscore_bed (with score in 5th col) provided, map using bedtools map
    if conscore_bed and bedtools_bin:
        # bedtools map -a tss_file -b conscore_bed -c 5 -o mean > outdir/breed_ch..._conscore.bed
        mapped = Path(outdir) / f"{breed}_chr6_conscore.bed"
        cmd = [bedtools_bin, "map", "-a", str(tss_file), "-b", conscore_bed, "-c", "5", "-o", "mean"]
        try:
            res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            print("[ERROR] bedtools map failed:", e)
            return None
        with open(mapped, 'w') as fh:
            fh.write(res.stdout)
        # cut fields 4,5 -> SNP, score
        with open(mapped) as fh, open(conscore_txt, 'w') as outfh:
            for ln in fh:
                parts = ln.rstrip('\n').split('\t')
                if len(parts) >= 5:
                    outfh.write(f"{parts[3]}\t{parts[4]}\n")
                elif len(parts) >= 4:
                    outfh.write(f"{parts[3]}\t0\n")
        print(f"[INFO] Created conscore: {conscore_txt}")
        return conscore_txt
    # else if conserved_bed + bedtools provided, do intersect -c (count) then produce conscore from count
    if conserved_bed and bedtools_bin:
        cons_bed_out = Path(outdir) / f"{breed}_chr6_conserved.bed"
        cmd = [bedtools_bin, "intersect", "-a", str(tss_file), "-b", conserved_bed, "-wa", "-c"]
        try:
            res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            print("[ERROR] bedtools intersect failed:", e)
            return None
        with open(cons_bed_out, 'w') as fh:
            fh.write(res.stdout)
        conscore_txt = Path(outdir) / f"{breed}_chr6_conscore.txt"
        with open(cons_bed_out) as fh, open(conscore_txt, 'w') as outfh:
            for ln in fh:
                parts = ln.rstrip('\n').split('\t')
                if len(parts) >= 5:
                    outfh.write(f"{parts[3]}\t{parts[4]}\n")
                elif len(parts) >= 4:
                    outfh.write(f"{parts[3]}\t0\n")
        print(f"[INFO] Created conscore (from intersect counts): {conscore_txt}")
        return conscore_txt
    print("[WARN] Could not produce conserved score; will fill Conserved/Conserved_score with 0")
    return None

def build_feature_file(prefix, outpath, breed, maf_path, l2_map, conscore_txt):
    # prefix.bim required
    bim = Path(prefix + ".bim")
    if not bim.exists():
        raise FileNotFoundError(f"BIM not found: {bim}")
    # read bim: columns chrom, snp, cm, pos, a1, a2 (no header)
    bim_df = pd.read_csv(bim, sep='\t', header=None, usecols=[0,1,4,5], names=['CHR','SNP','A1','A2'])
    # get MAF series
    maf_series = read_maf(maf_path, bim_df) if maf_path else pd.Series([0.0]*len(bim_df), index=bim_df['SNP'])
    # build l2 values (l2_map is dict variant->value, or None)
    if isinstance(l2_map, dict):
        l2_vals = [l2_map.get(snp, 0.0) for snp in bim_df['SNP']]
    else:
        l2_vals = [0.0]*len(bim_df)
    # conserved score
    cons_map = {}
    cons_score_map = {}
    if conscore_txt and Path(conscore_txt).exists():
        try:
            cons_df = pd.read_csv(conscore_txt, sep='\t', header=None, names=['SNP','Conserved_score'])
            cons_score_map = dict(zip(cons_df['SNP'].astype(str), cons_df['Conserved_score'].astype(float)))
            cons_map = {k:(1 if cons_score_map.get(k,0)>0 else 0) for k in cons_score_map.keys()}
        except Exception as e:
            print("[WARN] reading conscore failed:", e)
    # assemble
    out_df = bim_df.copy()
    out_df['MAF'] = [maf_series.get(snp, 0.0) for snp in out_df['SNP']]
    out_df['Conserved'] = [cons_map.get(snp, 0) for snp in out_df['SNP']]
    out_df['L2'] = l2_vals
    out_df['Conserved_score'] = [cons_score_map.get(snp, 0.0) for snp in out_df['SNP']]
    # enforce column order
    cols = ['CHR','SNP','A1','A2','MAF','Conserved','L2','Conserved_score']
    out_df = out_df[cols]
    out_df.to_csv(outpath, sep='\t', index=False)
    print(f"[INFO] Written feature file: {outpath}")
    return outpath

def read_maf(maf_path, bim_df):
    try:
        maf_df = pd.read_csv(maf_path, delim_whitespace=True, dtype={'SNP':str})
        maf_map = dict(zip(maf_df['SNP'].astype(str), maf_df['MAF'].astype(float)))
    except Exception as e:
        print(f"[WARN] failed to read maf {maf_path}: {e}")
        maf_map = {}
    return pd.Series([maf_map.get(snp, 0.0) for snp in bim_df['SNP']], index=bim_df['SNP'])

def main():
    p = argparse.ArgumentParser(description="Generate minimal feature file with MAF, Conserved, L2, Conserved_score")
    p.add_argument("--bfile", required=True, help="PLINK bfile prefix (no suffix), e.g. /path/TP_chr6")
    p.add_argument("--out", required=True, help="output feature file path")
    p.add_argument("--breed", required=True, help="breed name")
    p.add_argument("--conserved-bed", required=False, help="conserved regions bed (optional)")
    p.add_argument("--conscore-bed", required=False, help="conserved score bed (optional, 5th column score) ")
    p.add_argument("--plink", required=False, help="plink executable (optional, needed if maf missing)")
    p.add_argument("--ldsc", required=False, help="ldsc.py script path (optional, needed if l2 missing)")
    p.add_argument("--bedtools", required=False, default="bedtools", help="bedtools executable (default 'bedtools')")
    args = p.parse_args()

    prefix = args.bfile
    outpath = Path(args.out)
    outdir = outpath.parent
    outdir.mkdir(parents=True, exist_ok=True)
    breed = args.breed

    # sample name from bfile basename
    sample = Path(prefix).name

    # expected intermediate files in outdir (we prefer files in outdir to avoid polluting bfile dir)
    maf_expected = Path(outdir) / f"{sample}_maf.frq"
    l2_txt_expected = Path(outdir) / f"{sample}_l2.txt"
    ld_gz_expected = Path(outdir) / f"{sample}_ld.l2.ldscore.gz"
    conscore_expected = Path(outdir) / f"{breed}_chr6_conscore.txt"

    # 1) ensure maf
    maf_path = None
    # check in same dir as prefix first (common)
    candidate_maf1 = Path(prefix).with_suffix("")  # dummy
    # naming: some pipelines put maf under same folder as bfile with name {breed}_maf.frq
    candidate_maf = Path(prefix).parent / f"{sample}_maf.frq"
    if candidate_maf.exists():
        maf_path = candidate_maf
    elif maf_expected.exists():
        maf_path = maf_expected
    else:
        maf_path = ensure_maf(prefix, outdir, args.plink)

    # 2) ensure l2 map
    l2_map = None
    # check l2 txt in outdir or prefix dir
    if Path(prefix).parent.joinpath(f"{sample}_l2.txt").exists():
        df_l2 = pd.read_csv(Path(prefix).parent.joinpath(f"{sample}_l2.txt"), sep='\t', header=None, names=['SNP','L2'])
        l2_map = dict(zip(df_l2['SNP'].astype(str), df_l2['L2'].astype(float)))
    elif l2_txt_expected.exists():
        df_l2 = pd.read_csv(l2_txt_expected, sep='\t', header=None, names=['SNP','L2'])
        l2_map = dict(zip(df_l2['SNP'].astype(str), df_l2['L2'].astype(float)))
    elif ld_gz_expected.exists():
        # parse gz and create dict
        with gzip.open(ld_gz_expected, 'rt') as fh:
            fh.readline()
            l2map = {}
            for ln in fh:
                parts = ln.strip().split()
                if len(parts) >= 2:
                    var = parts[1]
                    try:
                        val = float(parts[-1])
                    except:
                        val = 0.0
                    l2map[var] = val
            l2_map = l2map
    else:
        l2_map = ensure_l2_txt(prefix, outdir, args.ldsc)

    # 3) ensure conserved score
    conscore = None
    # check candidate in outdir or prefix dir
    cand1 = Path(prefix).parent / f"{breed}_chr6_conscore.txt"
    if cand1.exists():
        conscore = cand1
    elif conscore_expected.exists():
        conscore = conscore_expected
    else:
        conscore = ensure_conserved_scores(prefix, outdir, breed, args.conserved_bed, args.conscore_bed, args.bedtools)

    # 4) build final feature
    build_feature_file(prefix, outpath, breed, maf_path, l2_map, conscore)

if __name__ == "__main__":
    main()
