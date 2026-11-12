import os
import sys
import argparse
import random
import traceback
import subprocess
import shutil
import hashlib
import multiprocessing
import glob
import time
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from sklearn.model_selection import KFold
import itertools
from datetime import datetime
import optuna
import json
import logging
import warnings
warnings.filterwarnings("ignore")
import threading
BASELINE_SEM = threading.Semaphore(1)
BASELINE_RESULTS = {}
BASELINE_THREADS = []
BASELINE_LOCK = threading.Lock()
LOG_LOCK = threading.Lock()
GLOBAL_TMP_DIR = None
def _baseline_runner(func, args, kwargs=None):
    with BASELINE_SEM:
        local_kwargs = dict(kwargs) if kwargs else {}
        result_key = local_kwargs.pop("_result_key", None)
        try:
            if local_kwargs:
                res = func(*args, **local_kwargs)
            else:
                res = func(*args)
        except Exception as e:
            if LOGGER:
                LOGGER.log(f"Baseline runner failed: {e}")
            res = None
        if result_key:
            with BASELINE_LOCK:
                BASELINE_RESULTS[result_key] = res
def launch_baseline(func, args, kwargs=None, log_name=None):
    t = threading.Thread(target=_baseline_runner, args=(func, args, kwargs))
    t.daemon = False
    t.start()
    BASELINE_THREADS.append(t)
    if LOGGER and log_name:
        LOGGER.log(f"Launching {log_name} in background")
    return t

# Configure logging
class Logger:
    def __init__(self, filename=None):
        if filename and os.path.exists(filename):
            try:
                if os.path.getsize(filename) > 0:
                    base, ext = os.path.splitext(filename)
                    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
                    archived = f"{base}.{ts}{ext if ext else '.log'}"
                    os.rename(filename, archived)
            except Exception:
                pass
        self.filename = filename
        if filename and os.path.dirname(filename):
            os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    def log(self, msg):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {msg}"
        with LOG_LOCK:
            print(line)
            if self.filename:
                with open(self.filename, "a") as f:
                    f.write(line + "\n")

LOGGER = None

def make_optuna_early_stop_callback(patience, min_delta):
    state = {"best": None, "bad": 0}
    patience = int(patience) if patience is not None else 0
    min_delta = float(min_delta) if min_delta is not None else 0.0

    def _cb(study, _trial):
        try:
            best = study.best_value
        except Exception:
            return
        if state["best"] is None or best > state["best"] + min_delta:
            state["best"] = best
            state["bad"] = 0
            return
        state["bad"] += 1
        if patience > 0 and state["bad"] >= patience:
            if LOGGER:
                LOGGER.log(f"Optuna early stop: no improvement for {patience} trials (min_delta={min_delta}).")
            try:
                study.stop()
            except Exception:
                return

    return _cb

try:
    import lightgbm as lgb
    HAS_LIGHTGBM = True
except Exception as e:
    lgb = None
    HAS_LIGHTGBM = False
    print(f"Warning: LightGBM import failed: {e}")

try:
    from pandas_plink import read_plink1_bin
    HAS_PANDAS_PLINK = True
except ImportError as e:
    HAS_PANDAS_PLINK = False
    print(f"Warning: pandas_plink import failed: {e}")
except Exception as e:
    HAS_PANDAS_PLINK = False
    print(f"Warning: pandas_plink import error: {e}")


from sklearn.metrics import mean_squared_error
from sklearn.linear_model import LinearRegression
from typing import Tuple
import multiprocessing
try:
    from threadpoolctl import threadpool_limits
except Exception:
    threadpool_limits = None

BLAS_THREADS = None
LGBM_N_JOBS = 1
PLINK_THREADS = 1
CLUMP_CACHE = {}
LD_PRUNE_CACHE = {}
CLUMP_MAX_SUMMARY_SNPS = None
INFO_EXTRACT_CACHE = {}
BIM_POS_CACHE = {}
SNP_SET_CACHE = {}
CLUMP_LD_CACHE = {}
CLUMP_CACHED_MAX_SNPS = None

def _read_snp_set(file_path: str):
    if not file_path:
        return None
    file_path = os.path.abspath(file_path)
    if file_path in SNP_SET_CACHE:
        return SNP_SET_CACHE[file_path]
    if not os.path.exists(file_path):
        SNP_SET_CACHE[file_path] = None
        return None
    try:
        s = set()
        with open(file_path, "r") as f:
            for line in f:
                tok = line.strip().split()
                if not tok:
                    continue
                s.add(tok[0])
        SNP_SET_CACHE[file_path] = s
        return s
    except Exception:
        SNP_SET_CACHE[file_path] = None
        return None

def _load_bim_positions(bfile: str):
    bfile = os.path.abspath(bfile)
    if bfile in BIM_POS_CACHE:
        return BIM_POS_CACHE[bfile]
    bim = f"{bfile}.bim"
    if not os.path.exists(bim):
        BIM_POS_CACHE[bfile] = {}
        return BIM_POS_CACHE[bfile]
    pos = {}
    try:
        df = pd.read_csv(bim, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
        for _, r in df.iterrows():
            snp = str(r["SNP"])
            try:
                bp = int(r["BP"])
            except Exception:
                continue
            chr_ = str(r["CHR"])
            pos[snp] = (chr_, bp)
    except Exception:
        pos = {}
    BIM_POS_CACHE[bfile] = pos
    return pos

def _parse_plink2_ld_file(ld_file: str, bfile: str):
    import gzip
    neighbors = {}
    pos_map = _load_bim_positions(bfile)
    opener = gzip.open if ld_file.endswith(".gz") else open
    with opener(ld_file, "rt") as f:
        header = f.readline().strip().split()
        if not header:
            return neighbors
        col = {c: i for i, c in enumerate(header)}
        def idx(*names):
            for n in names:
                if n in col:
                    return col[n]
            return None
        ia = idx("SNP_A", "ID_A", "ID1", "ID_A1", "SNP1")
        ib = idx("SNP_B", "ID_B", "ID2", "ID_B1", "SNP2")
        iba = idx("BP_A", "POS_A", "BP1", "POS1")
        ibb = idx("BP_B", "POS_B", "BP2", "POS2")
        ir2 = idx("R2", "R2_UNPHASED", "UNPHASED_R2", "PHASED_R2")
        if ia is None or ib is None or ir2 is None:
            return neighbors
        for line in f:
            parts = line.strip().split()
            if len(parts) <= max(ia, ib, ir2):
                continue
            a = parts[ia]
            b = parts[ib]
            if not a or not b or a == b:
                continue
            try:
                r2 = float(parts[ir2])
            except Exception:
                continue
            bpa = None
            bpb = None
            if iba is not None and iba < len(parts):
                try:
                    bpa = int(float(parts[iba]))
                except Exception:
                    bpa = None
            if ibb is not None and ibb < len(parts):
                try:
                    bpb = int(float(parts[ibb]))
                except Exception:
                    bpb = None
            if bpa is None or bpb is None:
                pa = pos_map.get(a)
                pb = pos_map.get(b)
                if pa is not None:
                    bpa = pa[1]
                if pb is not None:
                    bpb = pb[1]
            dist = abs(bpa - bpb) if (bpa is not None and bpb is not None) else 0
            neighbors.setdefault(a, []).append((b, r2, dist))
            neighbors.setdefault(b, []).append((a, r2, dist))
    return neighbors

def _build_clump_ld_cache_plink2(bfile: str, snps: list, out_prefix: str, max_kb: int, min_r2: float, plink_threads: int):
    prefix = os.path.abspath(out_prefix)
    snp_file = f"{prefix}.snps"
    ld_prefix = f"{prefix}.ld"
    ld_txt = f"{ld_prefix}.ld"
    ld_gz = f"{ld_prefix}.ld.gz"
    ld_vcor = f"{ld_prefix}.vcor"
    ld_vcor_gz = f"{ld_prefix}.vcor.gz"
    if not os.path.exists(snp_file):
        try:
            with open(snp_file, "w") as f:
                for s in snps:
                    f.write(str(s) + "\n")
        except Exception:
            return None
    if os.path.exists(ld_gz):
        return ld_gz
    if os.path.exists(ld_txt):
        return ld_txt
    if os.path.exists(ld_vcor_gz):
        return ld_vcor_gz
    if os.path.exists(ld_vcor):
        return ld_vcor
    cmd = (
        f"{PLINK2_EXEC} --bfile {bfile} --extract {snp_file} "
        f"--r2-unphased --ld-window-kb {int(max_kb)} --ld-window-r2 {float(min_r2)} --ld-window 999999 "
        f"--threads {int(plink_threads)} --allow-no-sex --out {ld_prefix}"
    )
    if not run_command(cmd):
        return None
    if os.path.exists(ld_vcor_gz):
        return ld_vcor_gz
    if os.path.exists(ld_vcor):
        return ld_vcor
    if os.path.exists(ld_gz):
        return ld_gz
    if os.path.exists(ld_txt):
        return ld_txt
    return None

def _get_or_build_clump_ld_cache(bfile: str, gwas_file: str, is_fastgwa: bool, out_prefix: str, p_max: float, max_kb: int, min_r2: float, plink_threads: int):
    top_n = int(CLUMP_CACHED_MAX_SNPS) if CLUMP_CACHED_MAX_SNPS else None
    base_key = (os.path.abspath(bfile), os.path.abspath(gwas_file), bool(is_fastgwa), top_n)
    cache = CLUMP_LD_CACHE.get(base_key)
    need = False
    if cache is None:
        need = True
    else:
        if float(p_max) > float(cache["p_max"]):
            need = True
        if int(max_kb) > int(cache["max_kb"]):
            need = True
        if float(min_r2) < float(cache["min_r2"]):
            need = True
    if not need:
        return cache
    p_max_eff = float(p_max) if cache is None else float(max(p_max, cache["p_max"]))
    max_kb_eff = int(max_kb) if cache is None else int(max(max_kb, cache["max_kb"]))
    min_r2_eff = float(min_r2) if cache is None else float(min(min_r2, cache["min_r2"]))

    pfile = ensure_two_col_gwas_file_filtered(gwas_file, is_fastgwa=is_fastgwa, p_max=p_max_eff, top_n=top_n)
    if not pfile or not os.path.exists(pfile):
        CLUMP_LD_CACHE[base_key] = None
        return None
    pairs = []
    try:
        with open(pfile, "r") as f:
            header = f.readline()
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    parts = line.rstrip("\n").split()
                    if len(parts) < 2:
                        continue
                snp = parts[0]
                try:
                    pv = float(parts[1])
                except Exception:
                    continue
                pairs.append((pv, snp))
    except Exception:
        return None
    if not pairs:
        return None
    pairs.sort(key=lambda x: x[0])
    pvals = {snp: pv for pv, snp in pairs}
    snps = [snp for _, snp in pairs]

    cache_key = (base_key, p_max_eff, max_kb_eff, min_r2_eff)
    cache_prefix = os.path.join(os.path.dirname(os.path.abspath(out_prefix)), f"ldcache_{_stable_key_hash(cache_key)}")
    ld_file = _build_clump_ld_cache_plink2(bfile=bfile, snps=snps, out_prefix=cache_prefix, max_kb=max_kb_eff, min_r2=min_r2_eff, plink_threads=plink_threads)
    if not ld_file:
        return None
    neighbors = _parse_plink2_ld_file(ld_file=ld_file, bfile=bfile)
    cache = {
        "p_max": p_max_eff,
        "max_kb": max_kb_eff,
        "min_r2": min_r2_eff,
        "pairs": pairs,
        "pvals": pvals,
        "neighbors": neighbors,
    }
    CLUMP_LD_CACHE[base_key] = cache
    return cache

def safe_accuracy(metrics: dict) -> float:
    try:
        v = metrics.get("accuracy", float("nan"))
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return 0.0
        return float(v)
    except Exception:
        return 0.0

def _stable_key_hash(key) -> str:
    import hashlib
    return hashlib.md5(repr(key).encode("utf-8")).hexdigest()[:12]
SUBSET_BFILE_CACHE = {}

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def _limit_blas(n_threads: int):
    if n_threads is None or n_threads <= 0 or threadpool_limits is None:
        class _Noop:
            def __enter__(self): return None
            def __exit__(self, exc_type, exc, tb): return False
        return _Noop()
    return threadpool_limits(limits=n_threads)

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
PLINK_EXEC = os.environ.get("IDEA_GP_PLINK", shutil.which("plink") or "plink")
PLINK2_EXEC = os.environ.get("IDEA_GP_PLINK2", shutil.which("plink2") or "plink2")
RSCRIPT_EXEC = os.environ.get("IDEA_GP_RSCRIPT", shutil.which("Rscript") or "Rscript")
HIBLUP_EXEC = os.environ.get("IDEA_GP_HIBLUP", shutil.which("hiblup") or "hiblup")

def run_command(cmd):
    """Run a shell command."""
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        if LOGGER:
            LOGGER.log(f"Error running command: {cmd}")
            LOGGER.log(f"STDOUT: {e.stdout.decode()}")
            LOGGER.log(f"STDERR: {e.stderr.decode()}")
        else:
            print(f"Error running command: {cmd}")
            print(f"STDOUT: {e.stdout.decode()}")
            print(f"STDERR: {e.stderr.decode()}")
        return False
    return True

def compute_hiblup_xrm_from_bfile(bfile: str, out_prefix: str, hiblup_exec: str = HIBLUP_EXEC, threads: int = 1) -> str:
    """
    Use HIBLUP to build an additive relationship matrix (XRM/GRM) from a PLINK bfile.
    Returns the XRM prefix to be used with --xrm (typically <prefix>.GA), or None on failure.
    """
    prefix = os.path.abspath(out_prefix)
    ga_prefix = f"{prefix}.GA"
    ga_bin = f"{ga_prefix}.bin"
    ga_id = f"{ga_prefix}.id"
    if os.path.exists(ga_bin) and os.path.exists(ga_id):
        return ga_prefix
    prune_in = compute_ld_pruned_snps(bfile, f"{prefix}_xrm_prune", plink_threads=PLINK_THREADS)
    if prune_in is None:
        if LOGGER: LOGGER.log("HIBLUP XRM: LD pruning failed.")
        return None
    pruned_bfile = f"{prefix}_xrm_pruned"
    cmd_prune = (
        f"{PLINK_EXEC} --bfile {bfile} --extract {prune_in} "
        f"--make-bed --out {pruned_bfile} --allow-no-sex"
    )
    if not run_command(cmd_prune):
        if LOGGER: LOGGER.log("HIBLUP XRM: failed to create pruned bfile.")
        return None
    cmd = (
        f"{hiblup_exec} --make-xrm --code-method 1 "
        f"--bfile {pruned_bfile} --add --thread {threads} --out {prefix}"
    )
    if not run_command(cmd):
        if LOGGER: LOGGER.log("HIBLUP XRM construction failed.")
        return None
    if not (os.path.exists(ga_bin) and os.path.exists(ga_id)):
        if LOGGER: LOGGER.log(f"HIBLUP XRM files not found at {ga_prefix}.bin / .id")
        return None
    return ga_prefix

def compute_ld_pruned_snps(bfile, out_prefix, window_kb=50, step=5, r2=0.2, plink_threads=1):
    """
    Compute LD-pruned SNP set for baseline using PLINK --indep-pairwise.
    Returns path to .prune.in file or None on failure.
    """
    indep_prefix = os.path.abspath(out_prefix)
    cmd = (
        f"{PLINK_EXEC} --bfile {bfile} "
        f"--indep-pairwise {window_kb} {step} {r2} "
        f"--allow-no-sex --threads {plink_threads} --out {indep_prefix}"
    )
    if not run_command(cmd):
        return None
    prune_in = f"{indep_prefix}.prune.in"
    if not os.path.exists(prune_in):
        return None
    return prune_in

def run_gwas_only(bfile, pheno_file, keep_ids_file, out_prefix, threads, use_fastgwa=True, background=False, reuse=False):
    """
    Run GWAS using PLINK2 (fast) or PLINK 1.9 (compatible).
    Returns: (output_file_path, is_plink2_bool, process_object_if_background)
    """
    global PLINK2_EXEC
 
    # Check for existing outputs first (reuse)
    existing_plink2 = glob.glob(f"{out_prefix}.*.glm.linear")
    existing_plink1 = os.path.exists(f"{out_prefix}.assoc.linear")
    if reuse and (existing_plink2 or existing_plink1):
        if existing_plink2:
            return existing_plink2[0], True, None
        return f"{out_prefix}.assoc.linear", False, None
    
    # Decide tool to run if not reusing
    run_plink2 = False
    if use_fastgwa and os.path.exists(PLINK2_EXEC):
        run_plink2 = True
        
    if run_plink2:
        # PLINK2
        cmd = [
            PLINK2_EXEC,
            "--bfile", bfile,
            "--keep", keep_ids_file,
            "--pheno", pheno_file,
            "--glm", "allow-no-covars", "omit-ref",
            "--threads", str(threads),
            "--out", out_prefix
        ]
        # PLINK2 output name is typically out_prefix.PhenoName.glm.linear
        # We return the pattern to match.
        out_file = f"{out_prefix}.*.glm.linear"
        tool_name = "PLINK2"
    else:
        # PLINK 1.9
        cmd = [
            PLINK_EXEC,
            "--bfile", bfile,
            "--keep", keep_ids_file,
            "--pheno", pheno_file,
            "--linear", "--allow-no-sex",
            "--threads", str(threads),
            "--out", out_prefix
        ]
        out_file = f"{out_prefix}.assoc.linear"
        tool_name = "PLINK1.9"

    if background:
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return out_file, run_plink2, p
        except Exception as e:
            if LOGGER: LOGGER.log(f"Failed to launch background {tool_name}: {e}")
            return None, False, None
    else:
        cmd_str = " ".join(cmd)
        if run_command(cmd_str):
            if run_plink2:
                files = glob.glob(out_file)
                if files:
                    return files[0], run_plink2, None
                else:
                    if LOGGER: LOGGER.log(f"PLINK2 output not found matching {out_file}")
                    return None, False, None
            return out_file, run_plink2, None
        else:
            return None, False, None

def run_gwas_and_clump(bfile, pheno_file, train_ids_file, out_prefix, p_thresh, clump_r2, clump_kb, plink_threads=1, extract_file=None, use_fastgwa=True):
    """
    1. Run GWAS (PLINK2 or PLINK1.9).
    2. Clump based on GWAS results.
    3. Extract SNPs.
    """
    # 1. Run GWAS
    gwas_out_prefix = f"{out_prefix}.gwas"
    gwas_result_file, is_plink2, _ = run_gwas_only(bfile, pheno_file, train_ids_file, gwas_out_prefix, plink_threads, use_fastgwa=use_fastgwa)
    
    if not gwas_result_file or not os.path.exists(gwas_result_file):
        if LOGGER: LOGGER.log(f"GWAS failed: {gwas_result_file} not found")
        return []

    # 2. Clump using the shared helper (which handles conversion if needed)
    return clump_from_gwas(bfile, gwas_result_file, out_prefix, p_thresh, clump_r2, clump_kb, plink_threads, extract_file, is_fastgwa=is_plink2)

def ensure_gwas(bfile, pheno_file, keep_ids_file, gwas_prefix, plink_threads=1, reuse=True):
    """
    Run GWAS once for a given train split and phenotype. Reuse if exists and reuse=True.
    Returns path to .assoc.linear file or None.
    """
    # Use run_gwas_only to support PLINK2 automatically
    # But ensure_gwas is typically used for baseline or legacy.
    # If we want to upgrade ensure_gwas to use PLINK2, we should.
    # But for now, let's keep it simple or redirect it.
    
    # If we want to reuse existing files, we need to check for both PLINK1.9 and PLINK2 outputs.
    # For simplicity, let's just use run_gwas_only logic.
    
    out_file, is_plink2, _ = run_gwas_only(bfile, pheno_file, keep_ids_file, gwas_prefix, plink_threads, use_fastgwa=True, background=False)
    return out_file

def clump_from_gwas(bfile, gwas_qassoc, out_prefix, p_thresh, clump_r2, clump_kb, plink_threads=1, extract_file=None, is_fastgwa=False):
    """
    Run PLINK clumping starting from an existing GWAS .qassoc file.
    If is_fastgwa=True (meaning PLINK2), it handles the .glm.linear format.
    Returns list of clumped lead SNP IDs.
    """
    key = (os.path.abspath(bfile), os.path.abspath(gwas_qassoc), float(p_thresh), clump_r2, clump_kb, os.path.abspath(extract_file) if extract_file else "None", is_fastgwa, int(CLUMP_MAX_SUMMARY_SNPS) if CLUMP_MAX_SUMMARY_SNPS else None)
    if key in CLUMP_CACHE:
        snps_cached = CLUMP_CACHE[key]
        return list(snps_cached) if snps_cached is not None else []
        
    clump_out = f"{out_prefix}.clump"
    
    target_qassoc = ensure_two_col_gwas_file_filtered(
        gwas_qassoc,
        is_fastgwa=is_fastgwa,
        p_max=p_thresh,
        top_n=CLUMP_MAX_SUMMARY_SNPS
    )
    if not target_qassoc:
        CLUMP_CACHE[key] = []
        return []

    cmd_clump = (f"{PLINK_EXEC} --bfile {bfile} --clump {target_qassoc} "
                 f"--clump-p1 {p_thresh} --clump-p2 {p_thresh} "
        f"--clump-r2 {clump_r2} --clump-kb {clump_kb} "
        f"--allow-no-sex --threads {plink_threads} --out {clump_out}")
    
    if extract_file:
        cmd_clump += f" --extract {extract_file}"

    clumped_file = f"{clump_out}.clumped"
    proc = subprocess.run(cmd_clump, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ok = proc.returncode == 0
    if not ok and not os.path.exists(clumped_file):
        if LOGGER:
            LOGGER.log(f"Error running command: {cmd_clump}")
            LOGGER.log(f"STDOUT: {proc.stdout.decode()}")
            LOGGER.log(f"STDERR: {proc.stderr.decode()}")
        else:
            print(f"Error running command: {cmd_clump}")
            print(f"STDOUT: {proc.stdout.decode()}")
            print(f"STDERR: {proc.stderr.decode()}")
        CLUMP_CACHE[key] = []
        return []
    if not os.path.exists(clumped_file):
        return []
    try:
        df = pd.read_csv(clumped_file, sep=r"\s+")
        if df.empty or "SNP" not in df.columns:
            CLUMP_CACHE[key] = []
            return []
        
        if "P" in df.columns:
            df = df.sort_values("P")
        snps = df["SNP"].tolist()
        CLUMP_CACHE[key] = list(snps)
        return snps
    except Exception as e:
        if LOGGER: LOGGER.log(f"Error reading clumped file: {e}")
        return []

def ensure_two_col_gwas_file(gwas_file: str, is_fastgwa: bool = False) -> str:
    if gwas_file is None:
        return None
    gwas_file = os.path.abspath(gwas_file)
    if gwas_file.endswith(".plink_fmt.temp"):
        return gwas_file if os.path.exists(gwas_file) else None
    key = (gwas_file, bool(is_fastgwa))
    out = f"{gwas_file}.{_stable_key_hash(key)}.plink_fmt.temp"
    if os.path.exists(out):
        return out
    if not os.path.exists(gwas_file):
        return None
    awk = (
        "awk 'BEGIN{OFS=\"\\t\"} "
        "NR==1{for(i=1;i<=NF;i++){if($i==\"SNP\"||$i==\"ID\")s=i; if($i==\"P\")p=i; if($i==\"TEST\")t=i;} "
        "print \"SNP\",\"P\"; next} "
        "{if(s>0 && p>0){if(t==0 || $t==\"ADD\"){if($p!=\"NA\" && $p!=\"\"){print $s,$p}}}}' "
    )
    cmd = f"{awk} {gwas_file} > {out}"
    if run_command(cmd):
        return out if os.path.exists(out) else None
    return None

def ensure_two_col_gwas_file_filtered(gwas_file: str, is_fastgwa: bool = False, p_max: float = None, top_n: int = None) -> str:
    if gwas_file is None:
        return None
    gwas_file = os.path.abspath(gwas_file)
    key = (gwas_file, bool(is_fastgwa), float(p_max) if p_max is not None else None, int(top_n) if top_n else None)
    out = f"{gwas_file}.{_stable_key_hash(key)}.plink_fmt.temp"
    if os.path.exists(out):
        return out
    if not os.path.exists(gwas_file):
        return None
    raw = f"{out}.raw"
    pmax_val = float(p_max) if p_max is not None else 0.0
    awk = (
        "awk -v pmax=" + str(pmax_val) + " 'BEGIN{OFS=\"\\t\"} "
        "NR==1{for(i=1;i<=NF;i++){if($i==\"SNP\"||$i==\"ID\")s=i; if($i==\"P\")p=i; if($i==\"TEST\")t=i;} "
        "print \"SNP\",\"P\"; next} "
        "{if(s>0 && p>0){if(t==0 || $t==\"ADD\"){if($p!=\"NA\" && $p!=\"\"){pv=$p+0; if(pmax<=0 || pv<=pmax){print $s,$p}}}}}' "
    )
    cmd = f"{awk} {gwas_file} > {raw}"
    if not run_command(cmd):
        return None
    if not os.path.exists(raw):
        return None
    if top_n and int(top_n) > 0:
        top_n = int(top_n)
        cmd2 = f"(head -n 1 {raw}; tail -n +2 {raw} | LC_ALL=C sort -gk2,2 | head -n {top_n}) > {out}"
        if not run_command(cmd2):
            try:
                os.replace(raw, out)
            except Exception:
                return raw
        else:
            try:
                os.remove(raw)
            except Exception:
                pass
        return out if os.path.exists(out) else None
    try:
        os.replace(raw, out)
    except Exception:
        return raw
    return out if os.path.exists(out) else None

def ld_prune_snps(bfile: str, keep_ids_file: str, out_prefix: str, ld_r2: float, window_kb: int, plink_threads: int = 1, extract_file: str = None) -> str:
    key = (
        os.path.abspath(bfile),
        os.path.abspath(keep_ids_file) if keep_ids_file else "None",
        float(ld_r2),
        int(window_kb),
        os.path.abspath(extract_file) if extract_file else "None",
    )
    if key in LD_PRUNE_CACHE:
        return LD_PRUNE_CACHE[key]
    prune_prefix = os.path.abspath(out_prefix)
    prune_in = f"{prune_prefix}.prune.in"
    if os.path.exists(prune_in):
        LD_PRUNE_CACHE[key] = prune_in
        return prune_in
    step = 50
    cmd = (
        f"{PLINK_EXEC} --bfile {bfile} --indep-pairwise {window_kb} {step} {ld_r2} "
        f"--allow-no-sex --threads {plink_threads} --out {prune_prefix}"
    )
    if keep_ids_file:
        cmd += f" --keep {keep_ids_file}"
    if extract_file:
        cmd += f" --extract {extract_file}"
    if not run_command(cmd):
        LD_PRUNE_CACHE[key] = None
        return None
    if not os.path.exists(prune_in):
        LD_PRUNE_CACHE[key] = None
        return None
    LD_PRUNE_CACHE[key] = prune_in
    return prune_in

def select_snps_from_gwas(
    bfile: str,
    gwas_file: str,
    keep_ids_file: str,
    out_prefix: str,
    p_thresh: float,
    ld_r2: float,
    window_kb: int,
    plink_threads: int = 1,
    extract_file: str = None,
    is_fastgwa: bool = False,
    selection_method: str = "clump",
    max_snps: int = None,
):
    if selection_method == "hybrid":
        top = clump_from_gwas(
            bfile, gwas_file, out_prefix, p_thresh, ld_r2, window_kb,
            plink_threads, extract_file=extract_file, is_fastgwa=is_fastgwa
        )
        if max_snps and len(top) >= max_snps:
            return top[:max_snps]
        bg_r2 = float(min(0.2, float(ld_r2))) if ld_r2 is not None else 0.2
        key = (
            os.path.abspath(bfile),
            os.path.abspath(keep_ids_file) if keep_ids_file else "None",
            float(bg_r2),
            int(window_kb),
            os.path.abspath(extract_file) if extract_file else "None",
        )
        prune_prefix = os.path.join(os.path.dirname(os.path.abspath(out_prefix)), f"ldprune_{_stable_key_hash(key)}")
        prune_in = ld_prune_snps(
            bfile=bfile,
            keep_ids_file=keep_ids_file,
            out_prefix=prune_prefix,
            ld_r2=bg_r2,
            window_kb=window_kb,
            plink_threads=plink_threads,
            extract_file=extract_file,
        )
        bg = []
        if prune_in and os.path.exists(prune_in):
            try:
                with open(prune_in, "r") as f:
                    bg = [line.strip() for line in f if line.strip()]
            except Exception:
                bg = []

        selected = []
        seen = set()
        for snp in top:
            if snp in seen:
                continue
            selected.append(snp)
            seen.add(snp)
        for snp in bg:
            if max_snps and len(selected) >= max_snps:
                break
            if snp in seen:
                continue
            selected.append(snp)
            seen.add(snp)
        return selected

    if selection_method == "clump_cached":
        cache = _get_or_build_clump_ld_cache(
            bfile=bfile,
            gwas_file=gwas_file,
            is_fastgwa=is_fastgwa,
            out_prefix=out_prefix,
            p_max=p_thresh,
            max_kb=window_kb,
            min_r2=ld_r2,
            plink_threads=plink_threads,
        )
        if not cache:
            return []
        extract_set = _read_snp_set(extract_file) if extract_file else None
        pvals = cache["pvals"]
        neighbors = cache["neighbors"]
        pairs = cache["pairs"]
        excluded = set()
        selected = []
        kb_bp = int(window_kb) * 1000
        for pv, snp in pairs:
            if pv > p_thresh:
                break
            if extract_set is not None and snp not in extract_set:
                continue
            if snp in excluded:
                continue
            selected.append(snp)
            if max_snps and len(selected) >= max_snps:
                break
            for nb, r2, dist in neighbors.get(snp, []):
                if dist > kb_bp:
                    continue
                if r2 < ld_r2:
                    continue
                if extract_set is not None and nb not in extract_set:
                    continue
                pnb = pvals.get(nb)
                if pnb is None or pnb > p_thresh:
                    continue
                excluded.add(nb)
        return selected
    if selection_method == "clump":
        snps = clump_from_gwas(
            bfile, gwas_file, out_prefix, p_thresh, ld_r2, window_kb,
            plink_threads, extract_file=extract_file, is_fastgwa=is_fastgwa
        )
        if max_snps and len(snps) > max_snps:
            return snps[:max_snps]
        return snps

    pfile = ensure_two_col_gwas_file(gwas_file, is_fastgwa=is_fastgwa)
    if not pfile:
        return []
    key = (
        os.path.abspath(bfile),
        os.path.abspath(keep_ids_file) if keep_ids_file else "None",
        float(ld_r2),
        int(window_kb),
        os.path.abspath(extract_file) if extract_file else "None",
    )
    prune_prefix = os.path.join(os.path.dirname(os.path.abspath(out_prefix)), f"ldprune_{_stable_key_hash(key)}")
    prune_in = ld_prune_snps(
        bfile=bfile,
        keep_ids_file=keep_ids_file,
        out_prefix=prune_prefix,
        ld_r2=ld_r2,
        window_kb=window_kb,
        plink_threads=plink_threads,
        extract_file=extract_file,
    )
    if not prune_in or not os.path.exists(prune_in):
        return []
    try:
        with open(prune_in, "r") as f:
            prune_set = set(line.strip() for line in f if line.strip())
    except Exception:
        return []
    if not prune_set:
        return []

    import heapq
    topk = max_snps if (max_snps and max_snps > 0) else None
    heap = []
    results = []
    try:
        with open(pfile, "r") as f:
            header = f.readline()
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    parts = line.rstrip("\n").split()
                    if len(parts) < 2:
                        continue
                snp = parts[0]
                if snp not in prune_set:
                    continue
                try:
                    p = float(parts[1])
                except Exception:
                    continue
                if p > p_thresh:
                    continue
                if topk:
                    heapq.heappush(heap, (-p, snp))
                    if len(heap) > topk:
                        heapq.heappop(heap)
                else:
                    results.append((p, snp))
    except Exception:
        return []

    if topk:
        pairs = [(-pp, snp) for pp, snp in heap]
        pairs.sort(key=lambda x: x[0])
        return [snp for _, snp in pairs]
    results.sort(key=lambda x: x[0])
    return [snp for _, snp in results]

def prepare_external_qassoc(external_assoc, assoc_fmt, out_prefix):
    """
    Convert external GWAS summary (e.g. GEMMA .assoc.txt) to PLINK .qassoc-like
    file with columns: CHR, SNP, BP, P.
    """
    qassoc = f"{out_prefix}.qassoc"
    if os.path.exists(qassoc):
        return qassoc
    try:
        df = pd.read_csv(external_assoc, sep=r"\s+")
        if assoc_fmt == "gemma":
            required = ["chr", "rs", "ps", "p_wald"]
            for c in required:
                if c not in df.columns:
                    raise ValueError(f"External GEMMA assoc missing column: {c}")
            
            data_dict = {
                "CHR": df["chr"],
                "SNP": df["rs"],
                "BP": df["ps"],
                "P": df["p_wald"]
            }
            # Optional columns for PRS
            if "beta" in df.columns: data_dict["BETA"] = df["beta"]
            if "allele1" in df.columns: data_dict["A1"] = df["allele1"] 
            
            out_df = pd.DataFrame(data_dict)
            
        elif assoc_fmt == "plink":
            required = ["CHR", "SNP", "BP", "P"]
            for c in required:
                if c not in df.columns:
                    raise ValueError(f"External PLINK assoc missing column: {c}")
            
            # Keep optional columns for PRS if they exist
            optional = ["A1", "BETA", "OR", "NMISS", "STAT"]
            keep_cols = required + [c for c in optional if c in df.columns]
            out_df = df[keep_cols].copy()
        else:
            raise ValueError(f"Unknown assoc_fmt: {assoc_fmt}")
        out_df.to_csv(qassoc, sep="\t", index=False)
        return qassoc
    except Exception as e:
        if LOGGER: LOGGER.log(f"Error preparing external qassoc: {e}")
        return None

def gwas_diagnostics(gwas_file: str, is_fastgwa: bool = False, max_chunks: int = None):
    thresholds = [1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2]
    counts = {t: 0 for t in thresholds}
    n = 0
    zeros = 0
    chi2_samples = []
    rng = np.random.default_rng(42)

    try:
        with open(gwas_file, "r") as f:
            header = (f.readline() or "").strip().split()
    except Exception:
        return None

    p_col = None
    for c in ["P", "p", "p_wald", "PVAL", "Pvalue"]:
        if c in header:
            p_col = c
            break
    if not p_col:
        return None

    t_col = None
    for c in ["T_STAT", "T", "STAT", "Z", "Z_STAT", "WALD_STAT"]:
        if c in header:
            t_col = c
            break

    usecols = [p_col] + ([t_col] if t_col else [])
    chunksize = 2_000_000
    read_chunks = 0
    try:
        file_size = os.path.getsize(gwas_file)
    except Exception:
        file_size = None

    if file_size and file_size >= 1_000_000_000:
        try:
            p_i = header.index(p_col) + 1
        except Exception:
            return None
        script = (
            "BEGIN{FS=\"[ \\t]+\"; OFS=\"\\t\"} "
            "NR==1{next} "
            "{"
            "pv=$p; if(pv==\"\"||pv==\"NA\"||pv==\"NaN\"||pv==\"nan\") next; n++; pv+=0; "
            "if(pv==0) z++; "
            "if(pv<1e-8) c1++; "
            "if(pv<1e-6) c2++; "
            "if(pv<1e-5) c3++; "
            "if(pv<1e-4) c4++; "
            "if(pv<1e-3) c5++; "
            "if(pv<1e-2) c6++; "
            "if(pv<5e-2) c7++; "
            "} "
            "END{print n, z, c1+0, c2+0, c3+0, c4+0, c5+0, c6+0, c7+0}"
        )
        try:
            out = subprocess.check_output(["awk", "-v", f"p={p_i}", script, gwas_file], text=True).strip()
            parts = out.split("\t")
            if len(parts) != 9:
                return None
            n = int(float(parts[0]))
            zeros = int(float(parts[1]))
            counts = {
                1e-8: int(float(parts[2])),
                1e-6: int(float(parts[3])),
                1e-5: int(float(parts[4])),
                1e-4: int(float(parts[5])),
                1e-3: int(float(parts[6])),
                1e-2: int(float(parts[7])),
                5e-2: int(float(parts[8])),
            }
        except Exception:
            return None

        return {
            "n": n,
            "p_col": p_col,
            "t_col": t_col,
            "p_zero": zeros,
            "counts": counts,
            "lambda_gc": None,
            "file": gwas_file,
            "is_fastgwa": bool(is_fastgwa),
        }

    try:
        for chunk in pd.read_csv(gwas_file, sep=r"\s+", usecols=usecols, chunksize=chunksize):
            p = chunk[p_col].to_numpy(dtype=float, copy=False)
            n += int(p.size)
            zeros += int((p == 0).sum())
            for t in thresholds:
                counts[t] += int((p < t).sum())
            if t_col:
                tstat = chunk[t_col].to_numpy(dtype=float, copy=False)
                k = min(2000, tstat.size)
                if k > 0:
                    idx = rng.choice(tstat.size, size=k, replace=False)
                    chi2 = (tstat[idx] ** 2)
                    chi2 = chi2[np.isfinite(chi2)]
                    if chi2.size:
                        chi2_samples.append(chi2)
            read_chunks += 1
            if max_chunks and read_chunks >= int(max_chunks):
                break
    except Exception:
        return None

    lambda_gc = None
    if chi2_samples:
        chi2_all = np.concatenate(chi2_samples)
        if chi2_all.size:
            med = float(np.median(chi2_all))
            lambda_gc = med / 0.4549364

    return {
        "n": n,
        "p_col": p_col,
        "t_col": t_col,
        "p_zero": zeros,
        "counts": counts,
        "lambda_gc": lambda_gc,
        "file": gwas_file,
        "is_fastgwa": bool(is_fastgwa),
    }

def infer_architecture_preset(diag: dict):
    if not diag or not diag.get("n"):
        return "hybrid"
    n = float(diag["n"])
    c = diag["counts"]
    n_p1e8 = float(c.get(1e-8, 0))
    n_p1e6 = float(c.get(1e-6, 0))
    n_p1e4 = float(c.get(1e-4, 0))
    lam = diag.get("lambda_gc", None)

    if (lam is not None and lam >= 1.25 and n_p1e4 >= 200000) or (n_p1e4 / n >= 5e-6 and n_p1e6 >= 50000):
        return "polygenic"
    if n_p1e8 <= 200 and n_p1e6 <= 2000 and (lam is None or lam <= 1.15):
        return "oligogenic"
    return "hybrid"

def biased_search_space(preset: str, fast_mode: bool, base_max_snps: int, base_clump_cap: int):
    preset = (preset or "auto").lower()
    if preset not in {"polygenic", "oligogenic", "hybrid"}:
        preset = "hybrid"

    if fast_mode:
        if preset == "polygenic":
            p_space = [1e-3, 0.01, 0.05, 0.2, 1.0]
            r2_space = [0.2, 0.5, 0.8]
            kb_space = [250, 500, 1000]
            max_snps = int(max(base_max_snps, 60000))
            clump_cap = int(max(base_clump_cap, 500000))
        elif preset == "oligogenic":
            p_space = [1e-8, 1e-6, 1e-5, 1e-4, 1e-3]
            r2_space = [0.05, 0.1, 0.2]
            kb_space = [500, 1000, 2000]
            max_snps = int(min(base_max_snps, 20000))
            clump_cap = int(max(base_clump_cap, 200000))
        else:
            p_space = [1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.2, 1.0]
            r2_space = [0.1, 0.2, 0.5]
            kb_space = [250, 500, 1000]
            max_snps = int(base_max_snps)
            clump_cap = int(max(base_clump_cap, 300000))
    else:
        if preset == "polygenic":
            p_space = [1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
            r2_space = [0.2, 0.5, 0.8]
            kb_space = [250, 500, 1000, 2000]
            max_snps = int(max(base_max_snps, 60000))
            clump_cap = int(max(base_clump_cap, 500000))
        elif preset == "oligogenic":
            p_space = [1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 0.01]
            r2_space = [0.05, 0.1, 0.2]
            kb_space = [500, 1000, 2000]
            max_snps = int(min(base_max_snps, 20000))
            clump_cap = int(max(base_clump_cap, 200000))
        else:
            p_space = [1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
            r2_space = [0.1, 0.2, 0.5, 0.8]
            kb_space = [250, 500, 1000, 2000]
            max_snps = int(base_max_snps)
            clump_cap = int(max(base_clump_cap, 300000))

    if max_snps < 1000:
        max_snps = 1000
    if clump_cap < 50000:
        clump_cap = 50000
    return {
        "preset": preset,
        "p_space": p_space,
        "r2_space": r2_space,
        "kb_space": kb_space,
        "max_snps": max_snps,
        "clump_cap": clump_cap,
    }

def load_prior(prior_arg: str):
    if not prior_arg:
        return {}
    try:
        if os.path.exists(prior_arg):
            with open(prior_arg, "r") as f:
                obj = json.load(f)
        else:
            obj = json.loads(prior_arg)
        return obj if isinstance(obj, dict) else {}
    except Exception:
        return {}

def infer_trait_name_from_pheno(pheno_file: str) -> str:
    if not pheno_file:
        return "trait"
    try:
        pheno_file = os.path.abspath(pheno_file)
    except Exception:
        pass
    if not os.path.exists(pheno_file):
        base = os.path.basename(pheno_file)
        stem = os.path.splitext(base)[0]
        return stem or "trait"
    try:
        with open(pheno_file, "r") as f:
            header = (f.readline() or "").strip().split()
        if header:
            if "FID" in header and "IID" in header:
                candidates = [c for c in header if c not in {"FID", "IID"}]
                if candidates:
                    return candidates[-1]
            return header[-1]
    except Exception:
        pass
    base = os.path.basename(pheno_file)
    stem = os.path.splitext(base)[0]
    return stem or "trait"

def apply_output_layout(args):
    if getattr(args, "no_organize_outputs", False):
        return
    raw_out = os.path.abspath(args.out)
    trait = infer_trait_name_from_pheno(args.pheno or getattr(args, "train_pheno", None))
    parts = [p for p in raw_out.split(os.sep) if p]
    while len(parts) >= 4 and parts[-4] == trait and parts[-2] == trait and parts[-3] == parts[-1]:
        parts = parts[:-2]
        raw_out = os.sep + os.sep.join(parts)
    run_root_dir = raw_out
    run_id = parts[-1] if parts else os.path.basename(raw_out)
    if len(parts) >= 3 and parts[-3] == trait and parts[-2] == parts[-1]:
        out_dir = os.sep + os.sep.join(parts[:-1])
        run_id = parts[-1]
    elif len(parts) >= 2 and parts[-2] == trait:
        out_dir = raw_out
        run_id = parts[-1]
    else:
        if len(parts) >= 2 and parts[-2] == parts[-1]:
            run_root_dir = os.path.dirname(raw_out)
            run_id = os.path.basename(raw_out)
        out_dir = os.path.join(run_root_dir, trait, run_id)
    try:
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(os.path.join(out_dir, "logs"), exist_ok=True)
        os.makedirs(os.path.join(out_dir, "folds"), exist_ok=True)
        os.makedirs(os.path.join(out_dir, "tmp"), exist_ok=True)
    except Exception:
        pass
    if not getattr(args, "log_file", None):
        args.log_file = os.path.join(out_dir, "logs", f"{run_id}.log")
    if not getattr(args, "tmp_dir", None):
        args.tmp_dir = os.path.join(out_dir, "tmp")
    args.out = os.path.join(out_dir, run_id)
    args.organized_out_dir = out_dir

def apply_prior_to_space(space: dict, prior: dict):
    if not space:
        return space
    if not prior:
        return space
    out = dict(space)
    if isinstance(prior.get("p_space"), list) and prior["p_space"]:
        out["p_space"] = [float(x) for x in prior["p_space"]]
    if isinstance(prior.get("r2_space"), list) and prior["r2_space"]:
        out["r2_space"] = [float(x) for x in prior["r2_space"]]
    if isinstance(prior.get("kb_space"), list) and prior["kb_space"]:
        out["kb_space"] = [int(x) for x in prior["kb_space"]]
    if prior.get("max_snps") is not None:
        try:
            out["max_snps"] = int(prior["max_snps"])
        except Exception:
            pass
    if prior.get("clump_cap") is not None:
        try:
            out["clump_cap"] = int(prior["clump_cap"])
        except Exception:
            pass
    if prior.get("preset") in {"polygenic", "oligogenic", "hybrid"}:
        out["preset"] = prior["preset"]
    return out

def choose_selection_method(preset: str, prior: dict = None):
    if prior and prior.get("selection_method") in {"clump", "prune", "hybrid", "clump_cached"}:
        return prior["selection_method"]
    preset = (preset or "hybrid").lower()
    if preset == "oligogenic":
        return "clump"
    if preset == "polygenic":
        return "prune"
    return "hybrid"

def train_lightgbm(train_bfile, test_bfile, train_pheno_file, test_pheno_file, train_ids_file, test_ids_file, snp_list_file):
    """
    Train LightGBM model. Handles separate train/test bfiles.
    """
    if not HAS_LIGHTGBM:
        if LOGGER: LOGGER.log("LightGBM not installed or import failed.")
        return 0.0

    temp_dir = GLOBAL_TMP_DIR if GLOBAL_TMP_DIR else os.path.dirname(train_ids_file)
    key = (os.path.abspath(train_bfile), os.path.abspath(test_bfile), os.path.abspath(snp_list_file))
    temp_prefix = os.path.join(temp_dir, f"temp_lgb_{abs(hash(key))}_{random.randint(1000,9999)}")
    
    # Helper to load data
    def load_data(bfile, pheno_file, ids_file, prefix_suffix):
        curr_prefix = f"{temp_prefix}_{prefix_suffix}"
        
        # Method 1: Pandas PLINK (Fast)
        if HAS_PANDAS_PLINK:
            try:
                if LOGGER: LOGGER.log(f"Loading {prefix_suffix} data using pandas_plink (fast)...")
                # Extract subset
                cmd_subset = (f"{PLINK_EXEC} --bfile {bfile} --extract {snp_list_file} "
                              f"--keep {ids_file} --make-bed --out {curr_prefix} "
                              f"--allow-no-sex --threads {PLINK_THREADS}")
                if not run_command(cmd_subset):
                    raise Exception("Subset failed")
                
                # Try reading with standard order (bim, fam, bed)
                try:
                    G = read_plink1_bin(f"{curr_prefix}.bim", f"{curr_prefix}.fam", f"{curr_prefix}.bed", verbose=False)
                except UnicodeDecodeError:
                    # Fallback for different versions: (bed, bim, fam)
                    G = read_plink1_bin(f"{curr_prefix}.bed", f"{curr_prefix}.bim", f"{curr_prefix}.fam", verbose=False)
                except Exception as e:
                     # If both fail, raise original error or let outer catch handle it
                     raise e
                
                X = G.transpose().values
                if hasattr(X, "compute"):
                    X = X.compute()
                
                fam = pd.read_csv(f"{curr_prefix}.fam", sep=r"\s+", header=None, usecols=[0, 1], names=["FID", "IID"])
                fam["FID"] = fam["FID"].astype(str)
                fam["IID"] = fam["IID"].astype(str)
                
                # Phenotype
                pheno_df = pd.read_csv(pheno_file, sep=r"\s+")
                if "Pheno" not in pheno_df.columns:
                    if isinstance(pheno_df.iloc[0, 2], str):
                        pheno_df = pd.read_csv(pheno_file, sep=r"\s+", header=0)
                    else:
                        pheno_df.columns = ["FID", "IID", "Pheno"]
                pheno_df["FID"] = pheno_df["FID"].astype(str)
                pheno_df["IID"] = pheno_df["IID"].astype(str)
                
                feat_names = [f"SNP_{i}" for i in range(X.shape[1])] # Note: this assumes same order if using same SNP list
                # Verify SNP order? PLINK --extract usually preserves order or sorts by pos. 
                # If we extract same list from different files, PLINK usually sorts by chr/pos.
                # To be safe, we should match by SNP ID.
                bim = pd.read_csv(f"{curr_prefix}.bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
                feat_names = bim["SNP"].tolist()
                
                df_X = pd.DataFrame(X, columns=feat_names)
                df_X["FID"] = fam["FID"]
                df_X["IID"] = fam["IID"]
                
                df_all = pd.merge(df_X, pheno_df, on=["FID", "IID"], how="inner")
                return df_all, feat_names
            except Exception as e:
                if LOGGER: LOGGER.log(f"Fast LightGBM load failed for {prefix_suffix}, falling back: {e}")
                pass

        # Method 2: Recode A (Slow)
        cmd_raw = (f"{PLINK_EXEC} --bfile {bfile} --extract {snp_list_file} "
                   f"--keep {ids_file} --recode A --out {curr_prefix} --allow-no-sex --threads {PLINK_THREADS}")
        if not run_command(cmd_raw):
            return None, None
        
        try:
            raw = pd.read_csv(f"{curr_prefix}.raw", sep=r"\s+")
            pheno_df = pd.read_csv(pheno_file, sep=r"\s+")
            if "Pheno" not in pheno_df.columns:
                if isinstance(pheno_df.iloc[0, 2], str):
                    pheno_df = pd.read_csv(pheno_file, sep=r"\s+", header=0)
                else:
                    pheno_df.columns = ["FID", "IID", "Pheno"]
            
            pheno_df["FID"] = pheno_df["FID"].astype(str)
            pheno_df["IID"] = pheno_df["IID"].astype(str)
            raw["FID"] = raw["FID"].astype(str)
            raw["IID"] = raw["IID"].astype(str)
            
            data = pd.merge(raw, pheno_df, on=["FID", "IID"])
            # raw columns: FID IID PAT MAT SEX PHENOTYPE SNP1_A ...
            feat_names = [c for c in raw.columns if c not in ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]]
            return data, feat_names
        except Exception as e:
            if LOGGER: LOGGER.log(f"Slow load failed: {e}")
            return None, None

    try:
        train_data, train_feats = load_data(train_bfile, train_pheno_file, train_ids_file, "train")
        test_data, test_feats = load_data(test_bfile, test_pheno_file, test_ids_file, "test")
        
        if train_data is None or test_data is None:
            return empty_metrics()
            
        # Align features
        common_feats = [f for f in train_feats if f in test_feats]
        if len(common_feats) == 0:
            if LOGGER: LOGGER.log("No common SNPs between train and test for LightGBM")
            return empty_metrics()
            
        X_train = train_data[common_feats]
        y_train = train_data["Pheno"]
        X_test = test_data[common_feats]
        y_test = test_data["Pheno"]
        
        params = {
            "objective": "regression",
            "metric": "rmse",
            "verbosity": -1,
            "boosting_type": "gbdt",
            "n_estimators": 300,
            "feature_fraction": 0.5,
            "bagging_fraction": 0.8,
            "bagging_freq": 1,
            "n_jobs": LGBM_N_JOBS,
        }
        model = lgb.LGBMRegressor(**params)
        model.fit(X_train, y_train)
        preds = model.predict(X_test)
        
        metrics = calculate_metrics(y_test, preds)
        if LOGGER: LOGGER.log(f"LightGBM Metrics: Acc={metrics['accuracy']:.4f}")
        return metrics
            
    except Exception as e:
        if LOGGER: LOGGER.log(f"LightGBM Error: {e}")
        return empty_metrics()
    finally:
        # Cleanup
        dir_name = os.path.dirname(temp_prefix)
        base = os.path.basename(temp_prefix)
        if os.path.exists(dir_name):
            for f in os.listdir(dir_name):
                if f.startswith(base):
                    try:
                        os.remove(os.path.join(dir_name, f))
                    except:
                        pass

def calculate_metrics(y_true, y_pred):
    """
    Calculate GS performance metrics.
    """
    if len(y_true) < 2:
        return {"accuracy": 0.0, "r2": 0.0, "mse": 0.0, "bias": 0.0, "spearman": 0.0}
    
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    
    # Pearson Correlation (Accuracy)
    acc = 0.0
    if np.std(y_true) > 0 and np.std(y_pred) > 0:
        acc = float(np.corrcoef(y_true, y_pred)[0, 1])
    
    # R2
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    
    # MSE
    mse = float(np.mean((y_true - y_pred) ** 2))
    
    # Bias (Slope of y_true ~ y_pred)
    bias = 0.0
    if np.var(y_pred) > 0:
        bias = np.cov(y_true, y_pred)[0, 1] / np.var(y_pred)
        
    # Spearman
    sp_rho = 0.0
    try:
        sp_rho, _ = spearmanr(y_true, y_pred)
        if np.isnan(sp_rho): sp_rho = 0.0
    except:
        pass
        
    return {
        "accuracy": acc,
        "r2": r2,
        "mse": mse,
        "bias": bias,
        "spearman": sp_rho
    }

def empty_metrics():
    return {"accuracy": float("nan"), "r2": float("nan"), "mse": float("nan"), "bias": float("nan"), "spearman": float("nan"), "h2": float("nan")}

def resolve_method(args_method: str, best_params: dict, proxy_metric: bool) -> str:
    if best_params and best_params.get("method"):
        return str(best_params["method"])
    if str(args_method) == "auto":
        return "rrblup"
    return str(args_method)

def ensure_info_extract(info_file: str, threshold: float, out_prefix: str, target_bfile: str = None) -> str:
    if not info_file or threshold is None:
        return None
    try:
        threshold = float(threshold)
    except Exception:
        return None
    if threshold <= 0:
        return None
    info_file = os.path.abspath(info_file)
    key = (info_file, threshold)
    cached = INFO_EXTRACT_CACHE.get(key)
    if cached and os.path.exists(cached):
        return cached
    out_dir = os.path.dirname(os.path.abspath(out_prefix))
    try:
        os.makedirs(out_dir, exist_ok=True)
    except Exception:
        pass
    out = os.path.join(out_dir, f"info_ge_{threshold:g}.{_stable_key_hash(key)}.snps")
    if os.path.exists(out):
        INFO_EXTRACT_CACHE[key] = out
        return out
    target_fmt = "unknown"
    if target_bfile:
        try:
            target_fmt = detect_id_format(target_bfile)
        except Exception:
            target_fmt = "unknown"
    if target_fmt not in {"chr_pos", "original"}:
        target_fmt = "unknown"

    awk = (
        "awk -v thr=" + str(threshold) + " '"
        "BEGIN{FS=\"[ \\t]+\"} "
        "NR==1{"
        "s=1;i=2;"
        "for(k=1;k<=NF;k++){"
        "c=$k; gsub(/\\r/,\"\",c);"
        "if(c==\"SNP\"||c==\"ID\"||c==\"rsid\"||c==\"RSID\")s=k;"
        "if(c==\"INFO\"||c==\"Info\"||c==\"DR2\"||c==\"R2\"||c==\"rsq\"||c==\"RSQ\")i=k;"
        "}"
        "next"
        "} "
        "{v=$i+0; if(v>=thr){id=$s; "
        + ("gsub(/^chr/,\"\",id); " if target_fmt == "original" else "")
        + ("if(id !~ /^chr/){id=\"chr\" id} " if target_fmt == "chr_pos" else "")
        + "if(!(id in seen)){seen[id]=1; print id}}}' "
    )
    cmd = f"{awk} {info_file} > {out}"
    if not run_command(cmd):
        return None
    if not os.path.exists(out):
        return None
    INFO_EXTRACT_CACHE[key] = out
    return out

def calculate_prs_proxy(train_bfile, test_bfile, train_ids_file, test_ids_file, test_pheno_file, snp_list_file, gwas_qassoc, out_prefix):
    """
    Calculate PRS as a fast proxy for GBLUP accuracy.
    Uses PLINK --score to calculate PRS on the validation set using GWAS effect sizes.
    Returns metric dictionary.
    """
    if not os.path.exists(gwas_qassoc):
        if LOGGER: LOGGER.log(f"PRS Proxy Error: GWAS file {gwas_qassoc} missing.")
        return calculate_metrics([], [])

    # Format of .qassoc: CHR SNP BP P OR/BETA ...
    # Format of .assoc.linear: CHR SNP BP A1 TEST NMISS BETA STAT P
    # Format of PLINK2 .glm.linear: #CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P
    
    score_file = f"{out_prefix}.score"
    try:
        # Read GWAS result to get Betas
        # Handle PLINK2 comment header if present
        df = pd.read_csv(gwas_qassoc, sep=r"\s+", comment=None)
        if df.columns[0].startswith("#"):
             df.rename(columns={df.columns[0]: df.columns[0][1:]}, inplace=True)
             
        # PLINK2 mappings
        if "ID" in df.columns and "SNP" not in df.columns:
            df.rename(columns={"ID": "SNP"}, inplace=True)
        
        # Filter for ADD test if present (standard output of --linear)
        if "TEST" in df.columns:
            df = df[df["TEST"] == "ADD"]
        
        # Identify Beta column
        beta_col = None
        for c in ["BETA", "OR", "log(OR)"]:
            if c in df.columns:
                beta_col = c
                break
        
        # Also need A1 (Effect Allele). .qassoc usually has 'A1'
        if "A1" not in df.columns or "SNP" not in df.columns or beta_col is None:
             if LOGGER: LOGGER.log(f"PRS Proxy Error: GWAS file missing A1, SNP, or BETA columns. Columns: {df.columns.tolist()}")
             return calculate_metrics([], [])
        
        # Filter by SNP list if provided
        if snp_list_file:
            target_snps = pd.read_csv(snp_list_file, header=None, names=["SNP"])
            df = df[df["SNP"].isin(target_snps["SNP"])]
            
        # Write score file: SNP A1 BETA
        df[["SNP", "A1", beta_col]].to_csv(score_file, sep="\t", index=False, header=False)
        
        # Run PLINK Score on Test Set
        # Note: --score <file> <snp-col> <allele-col> <score-col> header? No, we wrote no header.
        # Default: SNP(1), A1(2), Score(3) matches our write.
        
        prs_out = f"{out_prefix}.prs"
        cmd_score = (f"{PLINK_EXEC} --bfile {test_bfile} --score {score_file} "
                     f"--keep {test_ids_file} --out {prs_out} "
                     f"--allow-no-sex --threads {PLINK_THREADS}")
        
        if not run_command(cmd_score):
             return calculate_metrics([], [])
             
        # Read PRS Profile
        # FID IID PHENO CNT CNT2 SCORE
        if not os.path.exists(f"{prs_out}.profile"):
             return calculate_metrics([], [])
             
        prof = pd.read_csv(f"{prs_out}.profile", sep=r"\s+")
        
        # Load True Phenotypes
        # Note: .profile contains PHENO column, but it might be -9 if not in bfile.
        # Safer to load from test_pheno_file
        
        te_ph = pd.read_csv(test_pheno_file, sep=r"\s+")
        if "Pheno" not in te_ph.columns:
            if isinstance(te_ph.iloc[0, 2], str): te_ph = pd.read_csv(test_pheno_file, sep=r"\s+", header=0)
            else: te_ph.columns = ["FID", "IID", "Pheno"]
        
        te_ph["FID"] = te_ph["FID"].astype(str)
        te_ph["IID"] = te_ph["IID"].astype(str)
        prof["FID"] = prof["FID"].astype(str)
        prof["IID"] = prof["IID"].astype(str)
        
        merged = pd.merge(prof, te_ph, on=["FID", "IID"])
        
        y_true = merged["Pheno"].astype(float)
        y_pred = merged["SCORE"].astype(float) # The PRS
        
        metrics = calculate_metrics(y_true, y_pred)
        return metrics

    except Exception as e:
        if LOGGER: LOGGER.log(f"PRS Proxy Failed: {e}")
        return calculate_metrics([], [])

def calculate_sblup_proxy(train_bfile, test_bfile, train_ids_file, test_ids_file, test_pheno_file, snp_list_file, gwas_qassoc, out_prefix, gcta_path="gcta64"):
    """
    Calculate SBLUP (Summary-based BLUP) using GCTA as a proxy for GBLUP.
    Uses GWAS summary stats + Reference LD to estimate effects, then PLINK --score for prediction.
    """
    if not os.path.exists(gwas_qassoc):
        if LOGGER: LOGGER.log(f"SBLUP Proxy Error: GWAS file {gwas_qassoc} missing.")
        return calculate_metrics([], [])

    # 1. Prepare GCTA Summary Statistics Format (.ma)
    # Required: SNP A1 A2 freq b se p N
    freq_out = f"{out_prefix}.freq"
    cmd_freq = (f"{PLINK_EXEC} --bfile {train_bfile} --keep {train_ids_file} "
                f"--freq --allow-no-sex --threads {PLINK_THREADS} --out {freq_out}")
    if not run_command(cmd_freq):
         if LOGGER: LOGGER.log("SBLUP Proxy Error: Failed to calculate frequencies.")
         return calculate_metrics([], [])
         
    try:
        # Read Freq (SNP A1 A2 MAF NCHROBS)
        if os.path.exists(f"{freq_out}.frq"):
            df_freq = pd.read_csv(f"{freq_out}.frq", sep=r"\s+")
        else:
             if LOGGER: LOGGER.log("SBLUP Proxy Error: .frq file not found.")
             return calculate_metrics([], [])

        # Read GWAS
        df_gwas = pd.read_csv(gwas_qassoc, sep=r"\s+", comment=None)
        if df_gwas.columns[0].startswith("#"):
             df_gwas.rename(columns={df_gwas.columns[0]: df_gwas.columns[0][1:]}, inplace=True)
             
        # PLINK2 mappings
        if "ID" in df_gwas.columns and "SNP" not in df_gwas.columns:
            df_gwas.rename(columns={"ID": "SNP"}, inplace=True)
        if "OBS_CT" in df_gwas.columns and "NMISS" not in df_gwas.columns:
            df_gwas.rename(columns={"OBS_CT": "NMISS"}, inplace=True)
        
        # Merge
        # Rename GWAS columns to standard if needed
        # PLINK qassoc: SNP, BETA, SE, P, NMISS
        if "BETA" not in df_gwas.columns and "OR" in df_gwas.columns:
             df_gwas["BETA"] = np.log(df_gwas["OR"])
        
        df_merge = pd.merge(df_gwas, df_freq, on="SNP")
        
        # Columns for GCTA: SNP A1 A2 freq b se p N
        gcta_ma = pd.DataFrame({
            "SNP": df_merge["SNP"],
            "A1": df_merge["A1"],
            "A2": df_merge["A2"],
            "freq": df_merge["MAF"],
            "b": df_merge["BETA"],
            "se": df_merge["SE"],
            "p": df_merge["P"],
            "N": df_merge["NMISS"]
        })
        
        if snp_list_file:
             target_snps = pd.read_csv(snp_list_file, header=None, names=["SNP"])
             gcta_ma = gcta_ma[gcta_ma["SNP"].isin(target_snps["SNP"])]
        
        ma_file = f"{out_prefix}.ma"
        gcta_ma.to_csv(ma_file, sep="\t", index=False)
        
    except Exception as e:
        if LOGGER: LOGGER.log(f"SBLUP Proxy Error preparing .ma file: {e}")
        return calculate_metrics([], [])

    # 2. Run GCTA --cojo-sblup
    sblup_out_prefix = f"{out_prefix}.sblup_run"
    cmd_sblup = (f"{gcta_path} --bfile {train_bfile} --keep {train_ids_file} "
                 f"--cojo-sblup {ma_file} --cojo-wind 2000 "
                 f"--thread {PLINK_THREADS} --out {sblup_out_prefix}")
    
    if not run_command(cmd_sblup):
        # GCTA might fail if no SNPs passed
        return calculate_metrics([], [])
        
    # Output file: .sblup.cojo
    cojo_file = f"{sblup_out_prefix}.sblup.cojo"
    if not os.path.exists(cojo_file):
         return calculate_metrics([], [])
         
    # 3. Predict on Test Set using PLINK --score
    try:
        df_cojo = pd.read_csv(cojo_file, sep=r"\s+")
        score_file = f"{out_prefix}.sblup.score"
        df_cojo[["SNP", "A1", "b"]].to_csv(score_file, sep="\t", index=False, header=False)
        
        prs_out = f"{out_prefix}.sblup.pred"
        cmd_score = (f"{PLINK_EXEC} --bfile {test_bfile} --score {score_file} "
                     f"--keep {test_ids_file} --out {prs_out} "
                     f"--allow-no-sex --threads {PLINK_THREADS}")
                     
        if not run_command(cmd_score):
             return calculate_metrics([], [])
             
        # Read .profile
        profile_file = f"{prs_out}.profile"
        if not os.path.exists(profile_file):
             return calculate_metrics([], [])
             
        df_pred = pd.read_csv(profile_file, sep=r"\s+")
        
        # Read True Phenotypes
        df_pheno = pd.read_csv(test_pheno_file, sep=r"\s+")
        if "Pheno" not in df_pheno.columns:
             if isinstance(df_pheno.iloc[0, 2], str):
                  df_pheno = pd.read_csv(test_pheno_file, sep=r"\s+", header=0)
             else:
                  df_pheno.columns = ["FID", "IID", "Pheno"]
        
        df_pred["FID"] = df_pred["FID"].astype(str)
        df_pred["IID"] = df_pred["IID"].astype(str)
        df_pheno["FID"] = df_pheno["FID"].astype(str)
        df_pheno["IID"] = df_pheno["IID"].astype(str)
        
        merged = pd.merge(df_pheno, df_pred, on=["FID", "IID"])
        
        return calculate_metrics(merged["Pheno"], merged["SCORE"])
        
    except Exception as e:
        if LOGGER: LOGGER.log(f"SBLUP Proxy Error during prediction: {e}")
        return calculate_metrics([], [])

def train_hiblup_gblup(train_bfile, test_bfile, train_pheno_file, test_pheno_file, train_ids_file, test_ids_file, snp_list_file):
    """
    Train GBLUP using HIBLUP. Handles separate files by merging for GRM.
    Uses caching for GRM based on input hash.
    Returns dict of metrics.
    """
    temp_dir = os.path.dirname(train_ids_file)
    
    # --- Caching Logic ---
    def _hash_file_content(filepath):
        if not os.path.exists(filepath): return "missing"
        try:
            hasher = hashlib.md5()
            with open(filepath, 'rb') as f:
                hasher.update(f.read())
            return hasher.hexdigest()
        except:
            return "error"

    # Hash bfiles by path (assume constant), IDs by content (small, often temp), SNP list by path+mtime
    key_str = f"{os.path.abspath(train_bfile)}_{os.path.abspath(test_bfile)}"
    try:
        trp = os.path.abspath(train_pheno_file) if train_pheno_file else ""
        tep = os.path.abspath(test_pheno_file) if test_pheno_file else ""
        key_str += f"_{trp}"
        if trp and os.path.exists(trp):
            key_str += f"_{os.path.getmtime(trp)}"
        key_str += f"_{tep}"
        if tep and os.path.exists(tep):
            key_str += f"_{os.path.getmtime(tep)}"
    except Exception:
        pass
    
    # Optimize caching for CV: If the combined set of individuals is the same, reuse GRM.
    # Load IDs to create a canonical sorted string
    try:
        t1 = pd.read_csv(train_ids_file, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
        t2 = pd.read_csv(test_ids_file, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
        # Combine and sort
        all_ids = pd.concat([t1, t2]).drop_duplicates().sort_values(["FID", "IID"])
        # Create a hash of the combined IDs
        ids_str = all_ids.to_string(index=False, header=False)
        ids_hash = hashlib.md5(ids_str.encode()).hexdigest()
    except:
        ids_hash = f"{_hash_file_content(train_ids_file)}_{_hash_file_content(test_ids_file)}"

    key_str += f"_{ids_hash}"

    if snp_list_file:
        snp_abs = os.path.abspath(snp_list_file)
        key_str += f"_{snp_abs}"
        # Include mtime to detect content changes (important for random baselines)
        if os.path.exists(snp_abs):
            key_str += f"_{os.path.getmtime(snp_abs)}"
    else:
        key_str += "_all"
        
    file_hash = hashlib.md5(key_str.encode()).hexdigest()
    
    # Use a deterministic temp prefix based on hash
    # We can put it in a 'cache' folder to avoid clutter
    cache_dir = os.path.join(temp_dir, "hags_grm_cache")
    if not os.path.exists(cache_dir):
        try: os.makedirs(cache_dir, exist_ok=True)
        except: pass
        
    temp_prefix = os.path.join(cache_dir, f"grm_{file_hash}")
    metrics_cache_file = os.path.join(cache_dir, f"metrics_{file_hash}.json")
    if os.path.exists(metrics_cache_file):
        try:
            with open(metrics_cache_file, "r") as f:
                cached = json.load(f)
            if isinstance(cached, dict) and ("accuracy" in cached or "r2" in cached):
                if LOGGER:
                    LOGGER.log(f"Reusing cached GBLUP metrics: {metrics_cache_file}")
                return cached
        except Exception:
            pass
    
    # Check if GRM exists (HIBLUP outputs .xrm.bin and .xrm.id usually)
    grm_prefix = f"{temp_prefix}_merged"
    grm_exists = False
    
    # Check for common HIBLUP extensions
    if (os.path.exists(f"{grm_prefix}.xrm.bin") and os.path.exists(f"{grm_prefix}.xrm.id")):
        grm_exists = True
        grm_suffix = ".xrm"
    elif (os.path.exists(f"{grm_prefix}.GA.bin") and os.path.exists(f"{grm_prefix}.GA.id")):
        # Rename existing .GA files to .xrm for consistency
        try:
             shutil.move(f"{grm_prefix}.GA.bin", f"{grm_prefix}.xrm.bin")
             shutil.move(f"{grm_prefix}.GA.id", f"{grm_prefix}.xrm.id")
             grm_exists = True
             grm_suffix = ".xrm"
             if LOGGER: LOGGER.log(f"Renamed cached GRM: {grm_prefix}.GA -> .xrm")
        except Exception as e:
             if LOGGER: LOGGER.log(f"Failed to rename cached GRM: {e}")
             pass
    elif (os.path.exists(f"{grm_prefix}.xrm.ga") and os.path.exists(f"{grm_prefix}.xrm.id")):
        try:
             shutil.move(f"{grm_prefix}.xrm.ga", f"{grm_prefix}.xrm.bin")
             grm_exists = True
             grm_suffix = ".xrm"
        except: pass
    else:
        grm_suffix = ".xrm"
    
    if grm_exists and LOGGER:
         LOGGER.log(f"Reusing cached GRM: {grm_prefix} (suffix: {grm_suffix})")

    # Use a unique solve prefix to avoid conflicts
    solve_out = f"{temp_prefix}_solve_{random.randint(1000,9999)}"

    try:
        if not grm_exists:
            lock_path = f"{grm_prefix}.lock"

            def _acquire_grm_lock(timeout_sec=7200, stale_sec=21600, sleep_sec=2):
                start = time.time()
                while True:
                    try:
                        fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                        with os.fdopen(fd, "w") as f:
                            f.write(f"{os.getpid()} {time.time()}\n")
                        return True
                    except FileExistsError:
                        try:
                            if time.time() - os.path.getmtime(lock_path) > stale_sec:
                                os.remove(lock_path)
                                continue
                        except Exception:
                            pass
                        if time.time() - start > timeout_sec:
                            return False
                        time.sleep(sleep_sec)
                    except Exception:
                        if time.time() - start > timeout_sec:
                            return False
                        time.sleep(sleep_sec)

            got_lock = _acquire_grm_lock()
            if not got_lock:
                if LOGGER:
                    LOGGER.log(f"Timeout waiting for GRM lock: {lock_path}")
                return empty_metrics()

            try:
                if (os.path.exists(f"{grm_prefix}.xrm.bin") and os.path.exists(f"{grm_prefix}.xrm.id")):
                    grm_exists = True
                    grm_suffix = ".xrm"
                elif (os.path.exists(f"{grm_prefix}.GA.bin") and os.path.exists(f"{grm_prefix}.GA.id")):
                    try:
                        shutil.move(f"{grm_prefix}.GA.bin", f"{grm_prefix}.xrm.bin")
                        shutil.move(f"{grm_prefix}.GA.id", f"{grm_prefix}.xrm.id")
                        grm_exists = True
                        grm_suffix = ".xrm"
                        if LOGGER:
                            LOGGER.log(f"Renamed cached GRM: {grm_prefix}.GA -> .xrm")
                    except Exception:
                        pass

                if grm_exists:
                    if LOGGER:
                        LOGGER.log(f"Reusing cached GRM (post-lock): {grm_prefix} (suffix: {grm_suffix})")
                else:
                    same_bfile = (os.path.abspath(train_bfile) == os.path.abspath(test_bfile))
                    
                    if same_bfile:
                        all_ids_file = f"{temp_prefix}_all.ids"
                        with open(all_ids_file, 'w') as outfile:
                            with open(train_ids_file) as infile: outfile.write(infile.read())
                            if not train_ids_file.endswith('\n'): outfile.write('\n')
                            with open(test_ids_file) as infile: outfile.write(infile.read())
                        
                        cmd_merge = (f"{PLINK_EXEC} --bfile {train_bfile} " + 
                                     (f"--extract {snp_list_file} " if snp_list_file else "") + 
                                     f"--keep {all_ids_file} --make-bed --out {temp_prefix}_merged "
                                     f"--allow-no-sex --threads {PLINK_THREADS}")
                        
                        if not run_command(cmd_merge):
                             if LOGGER: LOGGER.log("PLINK extraction failed.")
                             return empty_metrics()
                    
                    else:
                        cmd_tr = (f"{PLINK_EXEC} --bfile {train_bfile} " + 
                                  (f"--extract {snp_list_file} " if snp_list_file else "") + 
                                  f"--keep {train_ids_file} --make-bed --out {temp_prefix}_train "
                                  f"--allow-no-sex --threads {PLINK_THREADS}")
                        if not run_command(cmd_tr): return empty_metrics()
                        
                        cmd_te = (f"{PLINK_EXEC} --bfile {test_bfile} " + 
                                  (f"--extract {snp_list_file} " if snp_list_file else "") + 
                                  f"--keep {test_ids_file} --make-bed --out {temp_prefix}_test "
                                  f"--allow-no-sex --threads {PLINK_THREADS}")
                        if not run_command(cmd_te): return empty_metrics()
                        
                        with open(f"{temp_prefix}_merge.txt", "w") as f:
                            f.write(f"{temp_prefix}_test.bed {temp_prefix}_test.bim {temp_prefix}_test.fam\n")
                            
                        cmd_merge = (f"{PLINK_EXEC} --bfile {temp_prefix}_train --merge-list {temp_prefix}_merge.txt "
                                     f"--make-bed --out {temp_prefix}_merged --allow-no-sex --threads {PLINK_THREADS}")
                        
                        if not run_command(cmd_merge):
                            if LOGGER: LOGGER.log("Merge for HIBLUP failed (likely ID conflict or SNP mismatch).")
                            return empty_metrics()
            
                    if os.path.exists(f"{temp_prefix}_merged.fam"):
                        with open(f"{temp_prefix}_merged.fam", 'r') as f:
                            num_inds = sum(1 for _ in f)
                        if LOGGER: LOGGER.log(f"HIBLUP Data Preparation: Merged dataset contains {num_inds} individuals.")
                        
                    if not (os.path.exists(f"{temp_prefix}_merged.bed") and os.path.exists(f"{temp_prefix}_merged.bim") and os.path.exists(f"{temp_prefix}_merged.fam")):
                        if LOGGER:
                            LOGGER.log(f"Error: merged PLINK files missing for GRM: {temp_prefix}_merged(.bed/.bim/.fam)")
                        return empty_metrics()

                    cmd_xrm = (f"{HIBLUP_EXEC} --make-xrm --bfile {temp_prefix}_merged "
                               f"--thread {BLAS_THREADS} --out {temp_prefix}_merged")
                    if not run_command(cmd_xrm):
                        for ext in [
                            ".GA.bin", ".GA.id",
                            ".xrm.bin", ".xrm.id",
                            ".log",
                        ]:
                            for p in [
                                f"{temp_prefix}_merged{ext}",
                                f"{temp_prefix}_merged.GA{ext}",
                            ]:
                                try:
                                    if os.path.exists(p):
                                        os.remove(p)
                                except Exception:
                                    pass
                        if not run_command(cmd_xrm):
                            return empty_metrics()
                    
                    if os.path.exists(f"{temp_prefix}_merged.GA.bin") and os.path.exists(f"{temp_prefix}_merged.GA.id"):
                         try:
                             shutil.move(f"{temp_prefix}_merged.GA.bin", f"{temp_prefix}_merged.xrm.bin")
                             shutil.move(f"{temp_prefix}_merged.GA.id", f"{temp_prefix}_merged.xrm.id")
                         except: pass
                         grm_suffix = ".xrm"
                    elif os.path.exists(f"{temp_prefix}_merged.xrm.bin") and os.path.exists(f"{temp_prefix}_merged.xrm.id"):
                         grm_suffix = ".xrm"
                    else:
                         grm_suffix = ".xrm"
            finally:
                try:
                    if os.path.exists(lock_path):
                        os.remove(lock_path)
                except Exception:
                    pass
        
        # 3. Prepare Phenotype (Train=Value, Test=NA)
        # Load Train Pheno
        tr_ph = pd.read_csv(train_pheno_file, sep=r"\s+")
        
        # Robust column detection
        if "FID" in tr_ph.columns and "IID" in tr_ph.columns:
            if "Pheno" in tr_ph.columns:
                target_col = "Pheno"
            else:
                # Use the last column that is not FID or IID
                cols = [c for c in tr_ph.columns if c not in ["FID", "IID"]]
                target_col = cols[-1] if cols else tr_ph.columns[-1]
            
            tr_ph = tr_ph[["FID", "IID", target_col]].copy()
            tr_ph.columns = ["FID", "IID", "Pheno"]
        else:
            # Fallback for files without standard headers
            if isinstance(tr_ph.iloc[0, 2], str): 
                # Retry with header=0 if it looks like we missed headers
                tr_ph = pd.read_csv(train_pheno_file, sep=r"\s+", header=0)
                if "FID" in tr_ph.columns and "IID" in tr_ph.columns:
                     if "Pheno" in tr_ph.columns: target_col = "Pheno"
                     else: target_col = tr_ph.columns[-1]
                     tr_ph = tr_ph[["FID", "IID", target_col]].copy()
                     tr_ph.columns = ["FID", "IID", "Pheno"]
                else:
                     # Still no standard headers, take 0, 1, and last
                     tr_ph = tr_ph.iloc[:, [0, 1, -1]]
                     tr_ph.columns = ["FID", "IID", "Pheno"]
            else:
                 # Assume 0, 1, and last column
                 tr_ph = tr_ph.iloc[:, [0, 1, -1]]
                 tr_ph.columns = ["FID", "IID", "Pheno"]
        
        tr_ids = pd.read_csv(train_ids_file, sep=r"\s+", header=None, names=["FID", "IID"])
        tr_ph["FID"] = tr_ph["FID"].astype(str); tr_ph["IID"] = tr_ph["IID"].astype(str)
        tr_ids["FID"] = tr_ids["FID"].astype(str); tr_ids["IID"] = tr_ids["IID"].astype(str)
        
        tr_data = tr_ids.merge(tr_ph, on=["FID", "IID"], how="inner")
        train_map = {k: v for k, v in zip(zip(tr_data["FID"], tr_data["IID"]), tr_data["Pheno"])}
        
        grm_id_file = f"{temp_prefix}_merged{grm_suffix}.id"
        if not os.path.exists(grm_id_file):
             if LOGGER: LOGGER.log(f"Error: GRM ID file missing: {grm_id_file}")
             return empty_metrics()
             
        # Detect ID file format (1 column vs 2 columns)
        grm_ids_df = pd.read_csv(grm_id_file, sep=r"\s+", header=None)
        
        hiblup_pheno = f"{solve_out}.pheno"
        
        if grm_ids_df.shape[1] == 1:
            # Single column ID (likely IID)
            grm_ids_df.columns = ["IID"]
            grm_ids_df["IID"] = grm_ids_df["IID"].astype(str)
            
            # Map IID -> Pheno
            # Warning: This assumes unique IIDs. If duplicates exist, last one wins or we need better logic.
            train_map_simple = {str(iid): pheno for iid, pheno in zip(tr_data["IID"], tr_data["Pheno"])}
            
            with open(hiblup_pheno, 'w') as f:
                # Write 3 columns: IID IID Pheno (Duplicate ID as FID to be safe for HIBLUP)
                # Using 2 columns (ID Pheno) might be risky if HIBLUP expects 3.
                # Safest bet for single-ID GRM is usually duplicating ID.
                # Use spaces, no header for HIBLUP if it's finicky, but usually it needs header.
                # Wait, "Error: non-digital phenotype value: 'YYTGXGG18100196'" suggests it's reading ID as phenotype.
                # This means it might be misinterpreting columns.
                # HIBLUP standard: FID IID Pheno.
                # If we provide FID IID Pheno, it should work.
                # Maybe the header is the issue? Or the separator?
                # Try NO HEADER and see.
                
                # No header
                for _, row in grm_ids_df.iterrows():
                    iid = str(row["IID"])
                    if iid in train_map_simple:
                        f.write(f"{iid} {train_map_simple[iid]}\n")
                    else:
                        f.write(f"{iid} NA\n")
                        
        else:
            # Two columns (FID IID)
            grm_ids_df.columns = ["FID", "IID"]
            grm_ids_df["FID"] = grm_ids_df["FID"].astype(str)
            grm_ids_df["IID"] = grm_ids_df["IID"].astype(str)
            
            # Map (FID, IID) -> Pheno
            train_map_full = {k: v for k, v in zip(zip(tr_data["FID"], tr_data["IID"]), tr_data["Pheno"])}
            
            with open(hiblup_pheno, 'w') as f:
                # No header for safety
                for _, row in grm_ids_df.iterrows():
                    key = (row["FID"], row["IID"])
                    if key in train_map_full:
                        f.write(f"{row['FID']} {row['IID']} {train_map_full[key]}\n")
                    else:
                        f.write(f"{row['FID']} {row['IID']} NA\n")

        
        # 4. Solve
        
        # Ensure .id and .bin files exist (HIBLUP v1.6.0+ might look for .id/.bin instead of .xrm.id/.xrm.bin)
        if os.path.exists(f"{temp_prefix}_merged.xrm.id") and not os.path.exists(f"{temp_prefix}_merged.id"):
             try: shutil.copy(f"{temp_prefix}_merged.xrm.id", f"{temp_prefix}_merged.id")
             except: pass
        if os.path.exists(f"{temp_prefix}_merged.xrm.bin") and not os.path.exists(f"{temp_prefix}_merged.bin"):
             try: shutil.copy(f"{temp_prefix}_merged.xrm.bin", f"{temp_prefix}_merged.bin")
             except: pass
             
        cmd_solve = (f"{HIBLUP_EXEC} --single-trait --pheno {hiblup_pheno} "
                     f"--xrm {temp_prefix}_merged --out {solve_out} "
                     f"--thread {BLAS_THREADS} --vc-method AI")
        if not run_command(cmd_solve): return empty_metrics()
        
        # 5. Extract Results and Validate
        rand_file = f"{solve_out}.rand"
        if not os.path.exists(rand_file): return empty_metrics()
        
        rand_df = pd.read_csv(rand_file, sep=r"\s+")
        
        # Normalize ID columns
        if "ID" in rand_df.columns:
            rand_df["ID"] = rand_df["ID"].astype(str)
            gen_col = None
            for c in rand_df.columns:
                if c != "ID" and not c.startswith("residual"):
                    gen_col = c
                    break
            if gen_col is None: return empty_metrics()
            
            # Load Test Truth
            te_ph = pd.read_csv(test_pheno_file, sep=r"\s+")
            
            # Robust column detection
            if "FID" in te_ph.columns and "IID" in te_ph.columns:
                 if "Pheno" in te_ph.columns: target_col = "Pheno"
                 else:
                     cols = [c for c in te_ph.columns if c not in ["FID", "IID"]]
                     target_col = cols[-1] if cols else te_ph.columns[-1]
                 te_ph = te_ph[["FID", "IID", target_col]].copy()
                 te_ph.columns = ["FID", "IID", "Pheno"]
            else:
                 if isinstance(te_ph.iloc[0, 2], str):
                      te_ph = pd.read_csv(test_pheno_file, sep=r"\s+", header=0)
                      if "FID" in te_ph.columns and "IID" in te_ph.columns:
                           if "Pheno" in te_ph.columns: target_col = "Pheno"
                           else: target_col = te_ph.columns[-1]
                           te_ph = te_ph[["FID", "IID", target_col]].copy()
                           te_ph.columns = ["FID", "IID", "Pheno"]
                      else:
                           te_ph = te_ph.iloc[:, [0, 1, -1]]
                           te_ph.columns = ["FID", "IID", "Pheno"]
                 else:
                      te_ph = te_ph.iloc[:, [0, 1, -1]]
                      te_ph.columns = ["FID", "IID", "Pheno"]
            te_ids = pd.read_csv(test_ids_file, sep=r"\s+", header=None, names=["FID", "IID"])
            te_ph["FID"] = te_ph["FID"].astype(str); te_ph["IID"] = te_ph["IID"].astype(str)
            te_ids["FID"] = te_ids["FID"].astype(str); te_ids["IID"] = te_ids["IID"].astype(str)
            
            test_truth = te_ids.merge(te_ph, on=["FID", "IID"], how="inner")
            
            # Try merging on IID -> ID
            final = test_truth.merge(rand_df, left_on="IID", right_on="ID", how="inner")
            
        elif "FID" in rand_df.columns and "IID" in rand_df.columns:
            rand_df["FID"] = rand_df["FID"].astype(str)
            rand_df["IID"] = rand_df["IID"].astype(str)
            gen_col = None
            for c in rand_df.columns:
                if c not in ["FID", "IID"] and not c.startswith("residual"):
                    gen_col = c
                    break
            if gen_col is None: return empty_metrics()
            
            # Load Test Truth
            te_ph = pd.read_csv(test_pheno_file, sep=r"\s+")
            
            # Robust column detection
            if "FID" in te_ph.columns and "IID" in te_ph.columns:
                 if "Pheno" in te_ph.columns: target_col = "Pheno"
                 else:
                     cols = [c for c in te_ph.columns if c not in ["FID", "IID"]]
                     target_col = cols[-1] if cols else te_ph.columns[-1]
                 te_ph = te_ph[["FID", "IID", target_col]].copy()
                 te_ph.columns = ["FID", "IID", "Pheno"]
            else:
                 if isinstance(te_ph.iloc[0, 2], str):
                      te_ph = pd.read_csv(test_pheno_file, sep=r"\s+", header=0)
                      if "FID" in te_ph.columns and "IID" in te_ph.columns:
                           if "Pheno" in te_ph.columns: target_col = "Pheno"
                           else: target_col = te_ph.columns[-1]
                           te_ph = te_ph[["FID", "IID", target_col]].copy()
                           te_ph.columns = ["FID", "IID", "Pheno"]
                      else:
                           te_ph = te_ph.iloc[:, [0, 1, -1]]
                           te_ph.columns = ["FID", "IID", "Pheno"]
                 else:
                      te_ph = te_ph.iloc[:, [0, 1, -1]]
                      te_ph.columns = ["FID", "IID", "Pheno"]
            te_ids = pd.read_csv(test_ids_file, sep=r"\s+", header=None, names=["FID", "IID"])
            te_ph["FID"] = te_ph["FID"].astype(str); te_ph["IID"] = te_ph["IID"].astype(str)
            te_ids["FID"] = te_ids["FID"].astype(str); te_ids["IID"] = te_ids["IID"].astype(str)
            
            test_truth = te_ids.merge(te_ph, on=["FID", "IID"], how="inner")
            
            final = test_truth.merge(rand_df, on=["FID", "IID"], how="inner")
            
        else:
            if LOGGER: LOGGER.log(f"HIBLUP Error: Unexpected columns in rand file: {rand_df.columns}")
            return empty_metrics()
        
        if len(final) < 2:
            if LOGGER: LOGGER.log(f"HIBLUP Warning: Only {len(final)} individuals matched in test set.")
            return empty_metrics()
        
        y_true = final["Pheno"].astype(float)
        y_pred = final[gen_col].astype(float)
        
        if y_true.std() == 0 or y_pred.std() == 0:
            if LOGGER: LOGGER.log("HIBLUP Warning: Zero variance in predictions or true values.")
            return empty_metrics()
        
        metrics = calculate_metrics(y_true, y_pred)
        
        # Log Variance Components and Calculate Heritability
        var_file = f"{solve_out}.var"
        if os.path.exists(var_file) and LOGGER:
            try:
                # Parse for Heritability
                var_df = pd.read_csv(var_file, sep=r"\s+")
                # Look for Variance components
                # Usually columns are: Component, Variance, SE
                # We need to find Va (Sigma2_a) and Ve (Sigma2_e)
                va = 0.0
                ve = 0.0
                
                # Check column names
                if "Variance" in var_df.columns:
                     # Filter by row names in first column
                     for idx, row in var_df.iterrows():
                         name = str(row.iloc[0])
                         val = float(row["Variance"])
                         if "Sigma2_a" in name or "Vg" in name:
                             va = val
                         elif "Sigma2_e" in name or "Ve" in name:
                             ve = val
                
                if va > 0 or ve > 0:
                    h2 = va / (va + ve) if (va + ve) > 0 else 0
                    metrics["h2"] = h2
                    LOGGER.log(f"  -> Heritability (h2): {h2:.4f}")
                    LOGGER.log(f"  -> Variance Components: Va={va:.4f}, Ve={ve:.4f}")
                else:
                    LOGGER.log("  -> Warning: Could not parse Va/Ve from HIBLUP output.")
                    if LOGGER: LOGGER.log(f"  -> Full Var File Content:\n{var_df.to_string()}")

            except Exception as e:
                LOGGER.log(f"Failed to parse heritability: {e}")
        
        if LOGGER: 
            LOGGER.log(f"  -> GBLUP Metrics: Acc={metrics['accuracy']:.4f}, R2={metrics['r2']:.4f}, Spear={metrics['spearman']:.4f}")
        try:
            if not os.path.exists(metrics_cache_file):
                tmp_path = f"{metrics_cache_file}.tmp.{os.getpid()}.{random.randint(1000,9999)}"
                with open(tmp_path, "w") as f:
                    json.dump(metrics, f, ensure_ascii=False)
                try:
                    os.replace(tmp_path, metrics_cache_file)
                except Exception:
                    try:
                        os.remove(tmp_path)
                    except Exception:
                        pass
        except Exception:
            pass
        return metrics
        
    except Exception as e:
        traceback.print_exc()
        if LOGGER: LOGGER.log(f"HIBLUP Trial Failed: {e}")
        return {"accuracy": 0.0, "r2": 0.0, "mse": 0.0, "bias": 0.0, "spearman": 0.0}
    finally:
        # Cleanup Solve files (but keep GRM cache)
        # solve_out files
        dir_name = os.path.dirname(solve_out) if 'solve_out' in locals() else None
        base = os.path.basename(solve_out) if 'solve_out' in locals() else None
        if dir_name and base and os.path.exists(dir_name):
            for f in os.listdir(dir_name):
                if f.startswith(base):
                    try:
                        os.remove(os.path.join(dir_name, f))
                    except:
                        pass
        # Cleanup temp PLINK files if not cached
        if not grm_exists and 'temp_prefix' in locals():
             # Usually we keep cached files.
             # But we should remove intermediate PLINK files (bed/bim/fam) if GRM is built.
             # HIBLUP needs .xrm.bin and .xrm.id, does it need original PLINK files for solve?
             # No, --xrm takes binary GRM.
             # So we can remove the temp_prefix PLINK files to save space in cache.
             for ext in [".bed", ".bim", ".fam", ".log", ".nosex", "_merge.txt"]:
                 # train/test/merged
                 for suffix in ["_train", "_test", "_merged"]:
                     f = f"{temp_prefix}{suffix}{ext}"
                     if os.path.exists(f):
                         try: os.remove(f)
                         except: pass
def evaluate_performance(train_bfile, test_bfile, train_pheno_file, test_pheno_file, train_ids_file, test_ids_file, snp_list_file, out_prefix, method="rrblup", libpath=None):
    if method in ["lightgbm", "adaptive", "hags"]:
        return train_lightgbm(train_bfile, test_bfile, train_pheno_file, test_pheno_file, train_ids_file, test_ids_file, snp_list_file)

    if method == "rrblup":
        exec_path = HIBLUP_EXEC if HIBLUP_EXEC else "hiblup"
        if exec_path and (os.path.isabs(exec_path) and os.path.exists(exec_path)):
            return train_hiblup_gblup(train_bfile, test_bfile, train_pheno_file, test_pheno_file, train_ids_file, test_ids_file, snp_list_file)
        if shutil.which(exec_path) is not None:
            return train_hiblup_gblup(train_bfile, test_bfile, train_pheno_file, test_pheno_file, train_ids_file, test_ids_file, snp_list_file)
        if LOGGER:
            LOGGER.log("HIBLUP not found; rrblup evaluation unavailable (R disabled).")
        return empty_metrics()

    if LOGGER:
        LOGGER.log(f"Unsupported method without R: {method}")
    return empty_metrics()

# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------

def load_info_scores(info_file):
    """
    Load imputation info scores (SNP, Info).
    Returns dict: {SNP: Info}
    """
    def _to_mapping(df: pd.DataFrame, snp_col: str, info_col: str):
        sub = df[[snp_col, info_col]].copy()
        sub[snp_col] = sub[snp_col].astype(str)
        sub[info_col] = pd.to_numeric(sub[info_col], errors="coerce")
        sub = sub.dropna(subset=[snp_col, info_col])
        return dict(zip(sub[snp_col], sub[info_col].astype(float)))

    try:
        df = pd.read_csv(info_file, sep=r"\s+", engine="python")
        cols = {str(c).strip().lower(): c for c in df.columns}

        snp_candidates = ["snp", "rs_id", "rsid", "id", "marker", "markername"]
        info_candidates = ["info", "impute_info", "rsq", "r2", "dr2"]

        snp_col = next((cols[k] for k in snp_candidates if k in cols), None)
        info_col = next((cols[k] for k in info_candidates if k in cols), None)
        if snp_col is not None and info_col is not None:
            return _to_mapping(df, snp_col, info_col)
    except Exception:
        pass

    try:
        df = pd.read_csv(info_file, sep=r"\s+", header=None, engine="python")
        if df.shape[1] < 2:
            return {}
        df = df.copy()
        df[0] = df[0].astype(str)
        numeric_cols = [c for c in df.columns[1:] if pd.to_numeric(df[c], errors="coerce").notna().any()]
        if not numeric_cols:
            return {}
        info_col = numeric_cols[-1]
        df[info_col] = pd.to_numeric(df[info_col], errors="coerce")
        df = df.dropna(subset=[0, info_col])
        return dict(zip(df[0], df[info_col].astype(float)))
    except Exception:
        return {}

def split_population(bfile, pheno_file, out_prefix, test_size=0.2):
    """
    Split population into Train and Test sets.
    Returns (train_bfile, test_bfile) paths.
    """
    from sklearn.model_selection import train_test_split
    
    # Read .fam file
    fam = pd.read_csv(f"{bfile}.fam", sep=r"\s+", header=None, usecols=[0,1], names=["FID", "IID"])
    # Read pheno to ensure we only split phenotyped individuals
    ph = pd.read_csv(pheno_file, sep=r"\s+")
    if "Pheno" not in ph.columns:
         if isinstance(ph.iloc[0, 2], str): ph = pd.read_csv(pheno_file, sep=r"\s+", header=0)
         else: ph.columns = ["FID", "IID", "Pheno"]
    
    # Merge
    merged = fam.merge(ph, on=["FID", "IID"])
    merged = merged[merged["Pheno"] != -9].dropna(subset=["Pheno"])
    
    train_ids, test_ids = train_test_split(merged[["FID", "IID"]], test_size=test_size, random_state=42)
    
    train_out = f"{out_prefix}.trial"
    test_out = f"{out_prefix}.test"
    
    train_ids.to_csv(f"{train_out}.ids", sep="\t", index=False, header=False)
    test_ids.to_csv(f"{test_out}.ids", sep="\t", index=False, header=False)
    
    # We don't create new bfiles here to save space/time, we just return the IDs file path?
    # No, the rest of the script expects bfile prefixes.
    # But wait, `evaluate_performance` takes `train_ids_file`.
    # So we can just return the original bfile + the ID files.
    # However, the script logic relies on `train_bfile` and `test_bfile`.
    # If we pass the SAME bfile but DIFFERENT IDs, it works (lines 865-866).
    
    return f"{train_out}.ids", f"{test_out}.ids"

def run_random_50k_baseline(train_bfile, test_bfile, train_ids, test_ids, train_pheno, test_pheno, out_prefix, info_file=None, info_thresh=0.3):
    """
    Run a quick baseline with random 50k SNPs (optionally filtered by Info).
    """
    # 1. Select SNPs
    if info_file:
        infos = load_info_scores(info_file)
        valid_snps = [s for s, i in infos.items() if i >= info_thresh]
    else:
        bim = pd.read_csv(f"{train_bfile}.bim", sep=r"\s+", header=None, usecols=[1], names=["SNP"])
        valid_snps = bim["SNP"].tolist()
        
    if len(valid_snps) > 50000:
        selected = random.sample(valid_snps, 50000)
    else:
        selected = valid_snps
        
    snp_file = f"{out_prefix}.random50k.snps"
    pd.DataFrame(selected, columns=["SNP"]).to_csv(snp_file, index=False, header=False)
    
    # 2. Run GBLUP
    metrics = train_hiblup_gblup(train_bfile, test_bfile, train_pheno, test_pheno, train_ids, test_ids, snp_file)
    return metrics

def run_fixed_snp_baseline(train_bfile, test_bfile, train_ids, test_ids, train_pheno, test_pheno, out_prefix, snp_list_file=None, description="Fixed Baseline"):
    """
    Run a baseline GBLUP with a fixed set of SNPs (or all if snp_list_file is None).
    """
    if LOGGER: LOGGER.log(f"Running Baseline: {description}")
    
    # Use train_hiblup_gblup directly
    # If snp_list_file is None, it runs on all SNPs in the bfile
    metrics = train_hiblup_gblup(train_bfile, test_bfile, train_pheno, test_pheno, train_ids, test_ids, snp_list_file)
    
    if LOGGER:
        LOGGER.log(f"Baseline {description} Results: Accuracy={metrics['accuracy']:.4f}, R2={metrics['r2']:.4f}")
        
    return metrics

def detect_id_format(bfile_prefix):
    """
    Detect if the bfile uses 'chr{CHR}_{POS}' format.
    """
    bim_file = f"{bfile_prefix}.bim"
    if not os.path.exists(bim_file):
        return "unknown"
    try:
        # Read first 10 lines
        df = pd.read_csv(bim_file, sep=r"\s+", header=None, nrows=10, usecols=[0, 1, 3], names=["CHR", "SNP", "BP"])
        
        is_chr_pos = True
        if len(df) == 0: return "unknown"
        
        for _, row in df.iterrows():
            # Check pattern chr{CHR}_{BP}
            expected = f"chr{row['CHR']}_{row['BP']}"
            if row["SNP"] != expected:
                is_chr_pos = False
                break
        return "chr_pos" if is_chr_pos else "original"
    except Exception as e:
        if LOGGER: LOGGER.log(f"Error detecting ID format for {bim_file}: {e}")
        return "unknown"

def get_snps_from_bim(bfile_prefix, target_format="original"):
    """
    Extract SNP IDs from a PLINK .bim file.
    If target_format is 'chr_pos', construct IDs as chr{CHR}_{BP}.
    """
    bim_file = f"{bfile_prefix}.bim"
    if not os.path.exists(bim_file):
        return []
    try:
        if target_format == "chr_pos":
             # Read CHR and BP
             # optimization: use dtypes
             df = pd.read_csv(bim_file, sep=r"\s+", header=None, usecols=[0, 3], names=["CHR", "BP"], engine="python")
             # Clean CHR column: remove 'chr' prefix if present to avoid duplication
             df["CHR"] = df["CHR"].astype(str).str.replace(r"^chr", "", case=False, regex=True)
             # Construct IDs: chr{CHR}_{BP}
             return ("chr" + df["CHR"] + "_" + df["BP"].astype(str)).tolist()
        else:
             # Read only 2nd column (SNP ID)
             # bim format: chr, snp, cm, bp, a1, a2
             df = pd.read_csv(bim_file, sep=r"\s+", header=None, usecols=[1], names=["SNP"], engine="python")
             return df["SNP"].tolist()
    except Exception as e:
        if LOGGER: LOGGER.log(f"Error reading bim file {bim_file}: {e}")
        return []

def read_snp_list_file(snp_list_file: str) -> list:
    if not snp_list_file or not os.path.exists(snp_list_file):
        return []
    snps = []
    try:
        with open(snp_list_file, "r") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                snps.append(s.split()[0])
    except Exception:
        return []
    return snps

def normalize_snp_list_for_bfile(snp_list_file: str, bfile_prefix: str, out_prefix: str) -> str:
    if not snp_list_file:
        return snp_list_file
    if not os.path.exists(snp_list_file):
        return snp_list_file
    fmt = detect_id_format(bfile_prefix)
    if fmt != "chr_pos":
        return snp_list_file

    needs = False
    try:
        with open(snp_list_file, "r") as f:
            for _ in range(50):
                line = f.readline()
                if not line:
                    break
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                tok = s.split()[0]
                if ":" in tok:
                    needs = True
                    break
    except Exception:
        return snp_list_file

    if not needs:
        return snp_list_file

    mapped = f"{out_prefix}.chrpos.snps"
    if os.path.exists(mapped):
        return mapped

    try:
        with open(snp_list_file, "r") as fin, open(mapped, "w") as fout:
            for line in fin:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                tok = s.split()[0]
                tok = tok.replace(":", "_")
                fout.write(tok + "\n")
        return mapped
    except Exception:
        return snp_list_file

def main():
    global PLINK_EXEC, PLINK2_EXEC, HIBLUP_EXEC, LOGGER, BLAS_THREADS, LGBM_N_JOBS, PLINK_THREADS
    parser = argparse.ArgumentParser(description="HAGS: Holo-Adaptive Genomic Selection Optimization")
    parser.add_argument("--bfile", default=None, help="PLINK binary prefix (if train/test are same)")
    parser.add_argument("--train_bfile", default=None, help="Training PLINK binary prefix")
    parser.add_argument("--test_bfile", default=None, help="Testing PLINK binary prefix")
    parser.add_argument("--pheno", default=None, help="Phenotype file (if train/test are same)")
    parser.add_argument("--train_pheno", default=None, help="Training Phenotype file")
    parser.add_argument("--test_pheno", default=None, help="Testing Phenotype file")
    
    parser.add_argument("--gwas_bfile", default=None, help="PLINK binary prefix for external GWAS (Population A)")
    parser.add_argument("--gwas_pheno", default=None, help="Phenotype file for external GWAS (Population A)")
    parser.add_argument("--out", required=True, help="Output prefix")
    parser.add_argument("--no_organize_outputs", action="store_true", default=False, help="Disable organized outputs")
    parser.add_argument("--folds", type=int, default=5, help="Number of CV folds")
    parser.add_argument("--fold_index", type=int, default=None, help="If set (1-based), only run this outer CV fold")
    parser.add_argument("--trials", type=int, default=20, help="Maximum number of optimization trials per fold")
    parser.add_argument("--method", type=str, default="auto", choices=["auto", "rrblup", "lasso", "kaml", "adaptive", "lightgbm"], help="GS method")
    parser.add_argument("--covar", type=str, default=None, help="Covariate file")
    parser.add_argument("--residualize", action="store_true", help="Residualize phenotype")
    parser.add_argument("--n_jobs", type=int, default=1, help="Parallel jobs for Optuna")
    parser.add_argument("--plink_threads", type=int, default=max(1, multiprocessing.cpu_count()//2), help="Threads for PLINK")
    parser.add_argument("--max_snps", type=int, default=30000, help="Cap selected SNPs")
    parser.add_argument("--reuse_gwas", action="store_true", help="Reuse GWAS", default=True)
    parser.add_argument("--blas_threads", type=int, default=max(1, multiprocessing.cpu_count()//2), help="Threads for BLAS")
    parser.add_argument("--lgbm_n_jobs", type=int, default=max(1, multiprocessing.cpu_count()//2), help="Threads for LightGBM")
    parser.add_argument("--external_assoc", type=str, default=None, help="External GWAS summary")
    parser.add_argument("--assoc_fmt", type=str, default="gemma", choices=["gemma", "plink"], help="Format of --external_assoc")
    parser.add_argument("--hiblup_exec", type=str, default=HIBLUP_EXEC, help="HIBLUP executable")
    parser.add_argument("--hiblup_xrm_prefix", type=str, default=None, help="HIBLUP XRM prefix")
    parser.add_argument("--baseline_bfile", type=str, default=None, help="PLINK bfile for HIBLUP baseline (single file)")
    parser.add_argument("--baseline_train_bfile", type=str, default=None, help="Baseline Training PLINK binary prefix")
    parser.add_argument("--baseline_test_bfile", type=str, default=None, help="Baseline Testing PLINK binary prefix")
    parser.add_argument("--hiblup_threads", type=int, default=max(1, multiprocessing.cpu_count()//2), help="Threads for HIBLUP")
    parser.add_argument("--optuna_patience", type=int, default=10, help="Patience")
    parser.add_argument("--optuna_min_delta", type=float, default=1e-3, help="Min delta")
    parser.add_argument("--optuna_seed", type=int, default=None, help="Seed")
    parser.add_argument("--optuna_startup_trials", type=int, default=10, help="Startup trials")
    parser.add_argument("--tmp_dir", type=str, default=None, help="Temp dir")
    parser.add_argument("--two_stage_opt", action="store_true", help="Two stage opt")
    parser.add_argument("--two_stage_topk", type=int, default=8, help="Top-K proxy trials to re-evaluate with full metric")
    parser.add_argument("--plink_exec", type=str, default=PLINK_EXEC, help="PLINK executable")
    parser.add_argument("--fast_mode", action="store_true", help="Fast mode")
    parser.add_argument("--skip_hiblup_baseline", action="store_true", help="Skip HIBLUP baseline")
    parser.add_argument("--lock_params_across_folds", action="store_true", help="Reuse best params")
    parser.add_argument("--exclude_ids", type=str, default=None, help="IDs to exclude")
    parser.add_argument("--libpath", type=str, default=None, help="Custom R lib path")
    parser.add_argument("--log_file", type=str, default=None, help="Log file path")
    parser.add_argument("--split_ratio", type=float, default=0.0, help="Ratio of test set (e.g. 0.2). If > 0, splits bfile into train/test.")
    parser.add_argument("--info_file", type=str, default=None, help="Imputation Info file (SNP Info)")
    parser.add_argument("--info_thresh", type=float, default=0.3, help="Default Info threshold")
    parser.add_argument("--proxy_metric", action="store_true", default=False, help="Use PRS/SBLUP as proxy metric in Inner CV instead of full GBLUP (faster)")
    parser.add_argument("--proxy_method", type=str, default="prs", choices=["prs", "sblup"], help="Method for proxy metric: 'prs' (C+T) or 'sblup' (GCTA-SBLUP)")
    parser.add_argument("--gcta_path", type=str, default=os.environ.get("IDEA_GP_GCTA", shutil.which("gcta64") or "gcta64"), help="Path to GCTA executable (required for --proxy_method sblup)")
    parser.add_argument("--inner_cv_folds", type=int, default=5, help="Number of folds for Inner CV (Trial)")
    parser.add_argument("--selection_method", type=str, default="auto", choices=["auto", "clump", "clump_cached", "prune", "hybrid"], help="SNP selection method")
    parser.add_argument("--optimize_selection_method", action="store_true", default=False, help="Tune SNP selection method (clump/prune/hybrid) in Optuna when --selection_method auto")
    parser.add_argument("--plink2_exec", type=str, default=PLINK2_EXEC, help="PLINK2 executable")
    parser.add_argument("--clump_max_summary_snps", type=int, default=200000, help="Cap rows passed to PLINK clump after p-filter")
    parser.add_argument("--clump_cache_max_snps", type=int, default=50000, help="Max SNPs used to build cached clump LD graph")
    parser.add_argument("--auto_search_space", action="store_true", default=False, help="Infer trait architecture from GWAS and bias search space")
    parser.add_argument("--architecture_preset", type=str, default="auto", choices=["auto", "polygenic", "oligogenic", "hybrid"], help="Manual override for trait architecture preset")
    parser.add_argument("--prior", type=str, default=None, help="JSON string or file path for prior settings (preset/selection_method/spaces/caps)")

    # New Baselines
    parser.add_argument("--baseline_full", action="store_true", help="Run full-genome GBLUP baseline")
    parser.add_argument("--baseline_custom_snps", type=str, default=None, help="SNP list for Custom Low-Density Baseline")
    parser.add_argument("--baseline_custom_bfile", type=str, default=None, help="PLINK bfile for Custom Low-Density Baseline (if different from input)")
    parser.add_argument("--baseline_commercial_snps", type=str, default=None, help="Comma-separated SNP lists for Commercial Chip Baselines")
    parser.add_argument("--baseline_commercial_bfile", type=str, default=None, help="PLINK bfile for Commercial Chip Baselines (if different from input)")
    parser.add_argument("--additional_baselines", type=str, default=None, help="JSON string or file path defining additional baselines. Format: [{'name': '...', 'bfile': '...', 'snp_list': '...'}, ...]")
    parser.add_argument("--baseline_only", action="store_true", default=False, help="Only compute baselines and update summary (skip Optuna/final evaluation)")
    parser.add_argument("--baseline_custom_plus_selected", action="store_true", default=False, help="Run a baseline on (custom chip SNPs ∪ selected SNPs) using input bfile")
    
    args = parser.parse_args()
    apply_output_layout(args)
    prior = load_prior(getattr(args, "prior", None))
    
    # Initialize Logger
    log_file = args.log_file if args.log_file else f"{args.out}.log"
    LOGGER = Logger(log_file)
    LOGGER.log(f"Starting HAGS Optimization with method: {args.method}")
    try:
        LOGGER.log(f"Resolved out_dir: {getattr(args, 'organized_out_dir', os.path.dirname(args.out))}")
        LOGGER.log(f"Resolved out_prefix: {args.out}")
    except Exception:
        pass
    LOGGER.log(f"Arguments: {args}")

    PLINK_EXEC = args.plink_exec
    PLINK2_EXEC = args.plink2_exec
    HIBLUP_EXEC = args.hiblup_exec
    BLAS_THREADS = args.blas_threads
    LGBM_N_JOBS = args.lgbm_n_jobs
    PLINK_THREADS = args.plink_threads
    global GLOBAL_TMP_DIR
    GLOBAL_TMP_DIR = args.tmp_dir
    global CLUMP_MAX_SUMMARY_SNPS
    CLUMP_MAX_SUMMARY_SNPS = args.clump_max_summary_snps
    global CLUMP_CACHED_MAX_SNPS
    CLUMP_CACHED_MAX_SNPS = args.clump_cache_max_snps
    if args.fast_mode and (CLUMP_MAX_SUMMARY_SNPS is None or CLUMP_MAX_SUMMARY_SNPS <= 0):
        CLUMP_MAX_SUMMARY_SNPS = 50000
    
    # Handle Train/Test files
    train_bfile = args.train_bfile if args.train_bfile else args.bfile
    test_bfile = args.test_bfile if args.test_bfile else (args.bfile if args.bfile else args.train_bfile)
    train_pheno = args.train_pheno if args.train_pheno else args.pheno
    test_pheno = args.test_pheno if args.test_pheno else (args.pheno if args.pheno else args.train_pheno)
    
    if not train_bfile or not train_pheno:
        LOGGER.log("Error: Missing genotype or phenotype file.")
        sys.exit(1)

    # Validate file existence
    if not os.path.exists(f"{train_bfile}.bed"):
        LOGGER.log(f"Error: Train bfile not found: {train_bfile}")
        sys.exit(1)
    if not os.path.exists(f"{test_bfile}.bed"):
        LOGGER.log(f"Error: Test bfile not found: {test_bfile}")
        sys.exit(1)

    global_info_scores = None
    if args.info_file:
        LOGGER.log(f"Info Scores enabled: {args.info_file}")

    # 1. Prepare Training/Testing IDs
    if args.split_ratio > 0 and args.bfile and not args.train_bfile and not args.test_bfile:
        LOGGER.log(f"Splitting population (ratio={args.split_ratio})...")
        trial_ids_path, test_ids_path = split_population(args.bfile, args.pheno, args.out, test_size=args.split_ratio)
        
        # Load these IDs into valid_train/valid_test
        valid_train = pd.read_csv(trial_ids_path, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
        valid_test = pd.read_csv(test_ids_path, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
        
        master_train_ids = trial_ids_path
        master_test_ids = test_ids_path
        
        LOGGER.log(f"Training Individuals: {len(valid_train)}")
        LOGGER.log(f"Testing Individuals: {len(valid_test)}")
        
    else:
        # Original Logic: Intersection of Bfile and Pheno
        # Prepare Training IDs
        train_fam = pd.read_csv(f"{train_bfile}.fam", sep=r"\s+", header=None, usecols=[0,1], names=["FID", "IID"])
        train_fam["FID"] = train_fam["FID"].astype(str)
        train_fam["IID"] = train_fam["IID"].astype(str)
        
        tr_ph = pd.read_csv(train_pheno, sep=r"\s+")
        # Check if we have enough columns for FID, IID, Pheno
        # Or if we should use specific columns (e.g. user wants specific trait)
        # For now, let's assume the input is corrected phenotypes with headers
        # FID, IID, Trait1, Trait2...
        # We need to select the trait to analyze. 
        # By default, let's take the LAST column if no specific trait is requested,
        # OR if there are 3 columns, take the 3rd.
        # But wait, the user provided a multi-column file.
        # Let's check headers.
        
        if "FID" in tr_ph.columns and "IID" in tr_ph.columns:
            # Good, we have IDs.
            # Now find the phenotype column.
            # If "Pheno" exists, use it.
            # Else, use the last column that is NOT FID or IID.
            if "Pheno" in tr_ph.columns:
                target_col = "Pheno"
            else:
                # Use the last column
                target_col = tr_ph.columns[-1]
                LOGGER.log(f"Auto-detected phenotype column: {target_col}")
            
            tr_ph = tr_ph[["FID", "IID", target_col]].copy()
            tr_ph.columns = ["FID", "IID", "Pheno"]
        else:
            # Fallback for no-header or simple format
            if tr_ph.shape[1] >= 3:
                # Assume first 3 are FID, IID, Pheno
                tr_ph = tr_ph.iloc[:, [0, 1, 2]]
                tr_ph.columns = ["FID", "IID", "Pheno"]
            else:
                 LOGGER.log("Error: Phenotype file must have at least 3 columns (FID, IID, Pheno)")
                 sys.exit(1)
        tr_ph["FID"] = tr_ph["FID"].astype(str)
        tr_ph["IID"] = tr_ph["IID"].astype(str)
        
        # Intersection
        valid_train = train_fam.merge(tr_ph, on=["FID", "IID"], how="inner")
        valid_train = valid_train[valid_train["Pheno"] != -9] # PLINK missing
        valid_train = valid_train.dropna(subset=["Pheno"])
        
        ids_dir = os.path.join(os.path.dirname(args.out), "ids")
        try:
            os.makedirs(ids_dir, exist_ok=True)
        except Exception:
            pass
        master_train_ids = os.path.join(ids_dir, "train.ids")
        valid_train[["FID", "IID"]].to_csv(master_train_ids, sep="\t", index=False, header=False)
        
        # Prepare Test IDs
        test_fam = pd.read_csv(f"{test_bfile}.fam", sep=r"\s+", header=None, usecols=[0,1], names=["FID", "IID"])
        test_fam["FID"] = test_fam["FID"].astype(str)
        test_fam["IID"] = test_fam["IID"].astype(str)
        
        te_ph = pd.read_csv(test_pheno, sep=r"\s+")
        
        if "FID" in te_ph.columns and "IID" in te_ph.columns:
            if "Pheno" in te_ph.columns:
                target_col = "Pheno"
            else:
                target_col = te_ph.columns[-1]
            
            te_ph = te_ph[["FID", "IID", target_col]].copy()
            te_ph.columns = ["FID", "IID", "Pheno"]
        else:
            if te_ph.shape[1] >= 3:
                te_ph = te_ph.iloc[:, [0, 1, 2]]
                te_ph.columns = ["FID", "IID", "Pheno"]
            else:
                 sys.exit(1)
        te_ph["FID"] = te_ph["FID"].astype(str)
        te_ph["IID"] = te_ph["IID"].astype(str)
        
        valid_test = test_fam.merge(te_ph, on=["FID", "IID"], how="inner")
        valid_test = valid_test[valid_test["Pheno"] != -9]
        valid_test = valid_test.dropna(subset=["Pheno"])
        
        master_test_ids = os.path.join(ids_dir, "test.ids")
        valid_test[["FID", "IID"]].to_csv(master_test_ids, sep="\t", index=False, header=False)
        
        LOGGER.log(f"Training Individuals: {len(valid_train)}")
        LOGGER.log(f"Testing Individuals: {len(valid_test)}")
    
    # If user provided external GWAS, prepare a reusable .qassoc once
    external_qassoc = None
    
    # New Logic: Generate External GWAS from provided GWAS bfile/pheno
    if args.gwas_bfile and args.gwas_pheno:
        LOGGER.log("Generating External GWAS summary from provided GWAS population...")
        if not os.path.exists(f"{args.gwas_bfile}.bed"):
            LOGGER.log(f"Error: GWAS bfile not found: {args.gwas_bfile}")
            sys.exit(1)
            
        gwas_out_prefix = f"{args.out}.external_gwas_run"
        # Create a temporary keep file for all individuals in GWAS bfile
        gwas_fam = pd.read_csv(f"{args.gwas_bfile}.fam", sep=r"\s+", header=None, usecols=[0,1], names=["FID", "IID"])
        gwas_keep = f"{gwas_out_prefix}.keep"
        gwas_fam.to_csv(gwas_keep, sep="\t", index=False, header=False)
        
        # Run GWAS
        qassoc_path = ensure_gwas(
            args.gwas_bfile, args.gwas_pheno, gwas_keep, gwas_out_prefix, 
            plink_threads=PLINK_THREADS, reuse=args.reuse_gwas
        )
        
        if qassoc_path and os.path.exists(qassoc_path):
            external_qassoc = qassoc_path
            LOGGER.log(f"Successfully generated External GWAS summary: {external_qassoc}")
        else:
            LOGGER.log("Failed to generate External GWAS summary. Aborting.")
            sys.exit(1)

    elif args.external_assoc:
        LOGGER.log(f"Using provided external GWAS summary: {args.external_assoc} (format={args.assoc_fmt})")
        external_qassoc = prepare_external_qassoc(args.external_assoc, args.assoc_fmt, f"{args.out}.external")
        if external_qassoc and os.path.exists(external_qassoc):
            LOGGER.log(f"Prepared external qassoc: {external_qassoc}")
        else:
            LOGGER.log("Failed to prepare external qassoc. Will fallback to internal GWAS.")
            external_qassoc = None
    
    # Parse Additional Baselines
    additional_baselines_list = []
    if args.additional_baselines:
        try:
            if args.additional_baselines.strip().startswith("[") or args.additional_baselines.endswith(".json"):
                # JSON format
                if args.additional_baselines.endswith(".json"):
                    with open(args.additional_baselines, 'r') as f:
                        additional_baselines_list = json.load(f)
                else:
                    additional_baselines_list = json.loads(args.additional_baselines)
            else:
                # String format: Name|Bfile|SNPList;...
                entries = args.additional_baselines.split(";")
                for entry in entries:
                    parts = entry.split("|")
                    if len(parts) >= 1:
                        name = parts[0].strip()
                        bfile = parts[1].strip() if len(parts) > 1 and parts[1].strip() else None
                        snp_list = parts[2].strip() if len(parts) > 2 and parts[2].strip() else None
                        additional_baselines_list.append({
                            "name": name,
                            "bfile": bfile,
                            "snp_list": snp_list
                        })
            LOGGER.log(f"Loaded {len(additional_baselines_list)} additional baselines.")
        except Exception as e:
            LOGGER.log(f"Error parsing --additional_baselines: {e}")

    # GS internal CV mode: if train and test use the same bfile, perform outer K-fold CV within GS
    if os.path.abspath(train_bfile) == os.path.abspath(test_bfile) and args.split_ratio == 0.0:
        LOGGER.log("Detected GS Internal CV mode (train/test bfile identical).")
        LOGGER.log(f"Strategy: {args.folds}-Fold CV within GS Population.")
        if external_qassoc:
             LOGGER.log("Feature Selection: Using External GWAS results (Population A).")
        else:
            LOGGER.log("Feature Selection: Preparing GWAS for Training Fold（优先复用已存在结果，缺失的将后台补算）。")

        run_dir = os.path.dirname(args.out)
        run_id = os.path.basename(args.out)
        summary_dir = os.path.join(run_dir, "summary")
        summary_txt = os.path.join(summary_dir, "cv_summary.txt")
        summary_json = os.path.join(summary_dir, "cv_summary.json")
        fold_records = []
        try:
            if os.path.exists(summary_json):
                with open(summary_json, "r") as f:
                    existing_summary = json.load(f)
                if isinstance(existing_summary, dict):
                    existing_folds = list(existing_summary.get("folds") or [])
                    done = []
                    missing = []
                    missing_baseline = []
                    for i in range(1, int(existing_summary.get("outer_folds", args.folds)) + 1):
                        fr = None
                        for x in existing_folds:
                            if isinstance(x, dict) and int(x.get("outer_fold", -1)) == int(i):
                                fr = x
                                break
                        if isinstance(fr, dict) and isinstance(fr.get("final_acc"), (int, float)) and not (isinstance(fr.get("final_acc"), float) and np.isnan(fr.get("final_acc"))):
                            done.append(i)
                        else:
                            missing.append(i)
                        if isinstance(fr, dict):
                            need_keys = ["baseline_full_acc", "baseline_custom_acc", "baseline_custom_plus_selected_acc", "baseline_commercial_acc"]
                            if any(fr.get(k) is None for k in need_keys):
                                missing_baseline.append(i)
                    if done or missing:
                        LOGGER.log(f"Existing summary status: completed_folds={len(done)}/{int(existing_summary.get('outer_folds', args.folds))}; missing_folds={missing}")
                    if missing_baseline:
                        LOGGER.log(f"Existing summary status: folds_missing_baselines={sorted(set(missing_baseline))}")
        except Exception:
            pass

        if args.baseline_only:
            try:
                with BASELINE_LOCK:
                    BASELINE_RESULTS.clear()
                del BASELINE_THREADS[:]
            except Exception:
                pass

            outer_accuracies = []
            mean_outer = float("nan")
            folds_n = int(args.folds)

            existing = None
            if os.path.exists(summary_json):
                try:
                    with open(summary_json, "r") as f:
                        existing = json.load(f)
                except Exception:
                    existing = None

            if isinstance(existing, dict):
                try:
                    folds_n = int(existing.get("outer_folds", folds_n))
                except Exception:
                    pass
                try:
                    outer_accuracies = list(existing.get("accuracies") or [])
                except Exception:
                    outer_accuracies = []
                # Guard against historical summaries containing null/non-numeric entries.
                outer_accuracies = [
                    float(x) for x in outer_accuracies
                    if isinstance(x, (int, float)) and not (isinstance(x, float) and np.isnan(x))
                ]
                if outer_accuracies:
                    try:
                        mean_outer = float(existing.get("mean_accuracy", float(np.mean(outer_accuracies))))
                    except Exception:
                        mean_outer = float(np.mean(outer_accuracies))
                try:
                    fold_records = list(existing.get("folds") or [])
                except Exception:
                    fold_records = []

            if not fold_records:
                fold_records = [{"outer_fold": int(i)} for i in range(1, folds_n + 1)]

            fold_by_no = {}
            for fr in fold_records:
                try:
                    fold_by_no[int(fr.get("outer_fold"))] = fr
                except Exception:
                    continue

            outer_kf = KFold(n_splits=folds_n, shuffle=True, random_state=42)
            X_indices = np.arange(len(valid_train))
            outer_splits = list(outer_kf.split(X_indices))

            for outer_i in range(1, folds_n + 1):
                if args.fold_index and outer_i != args.fold_index:
                    continue

                outer_dir = os.path.join(run_dir, "folds", f"outer{outer_i}")
                try:
                    os.makedirs(outer_dir, exist_ok=True)
                    os.makedirs(os.path.join(outer_dir, "ids"), exist_ok=True)
                except Exception:
                    pass

                ot_ids = os.path.join(outer_dir, "ids", "train.ids")
                ov_ids = os.path.join(outer_dir, "ids", "test.ids")

                if not (os.path.exists(ot_ids) and os.path.exists(ov_ids)):
                    try:
                        tr_idx, te_idx = outer_splits[outer_i - 1]
                        outer_train = valid_train.iloc[tr_idx]
                        outer_test = valid_train.iloc[te_idx]
                        outer_train[["FID", "IID"]].to_csv(ot_ids, sep="\t", index=False, header=False)
                        outer_test[["FID", "IID"]].to_csv(ov_ids, sep="\t", index=False, header=False)
                    except Exception as e:
                        LOGGER.log(f"Baseline-only: failed to prepare fold IDs for outer{outer_i}: {e}")
                        continue

                fr = fold_by_no.get(int(outer_i))
                if fr is None:
                    fr = {"outer_fold": int(outer_i)}
                    fold_records.append(fr)
                    fold_by_no[int(outer_i)] = fr

                try:
                    if "train_n" not in fr:
                        with open(ot_ids, "r") as f:
                            fr["train_n"] = sum(1 for _ in f)
                    if "test_n" not in fr:
                        with open(ov_ids, "r") as f:
                            fr["test_n"] = sum(1 for _ in f)
                except Exception:
                    pass

                baseline_bfile_arg = args.baseline_train_bfile if args.baseline_train_bfile else args.baseline_bfile
                baseline_dir = os.path.join(outer_dir, "baselines")
                try:
                    os.makedirs(baseline_dir, exist_ok=True)
                except Exception:
                    pass

                if not args.skip_hiblup_baseline and baseline_bfile_arg and fr.get("baseline_orig_acc") is None:
                    try:
                        target_id_fmt = detect_id_format(train_bfile)
                        extracted_snps = get_snps_from_bim(baseline_bfile_arg, target_format=target_id_fmt)
                        if extracted_snps:
                            temp_snp_file = os.path.join(baseline_dir, "baseline_orig.snps")
                            pd.DataFrame(extracted_snps, columns=["SNP"]).to_csv(temp_snp_file, index=False, header=False)
                            try:
                                launch_baseline(
                                    run_fixed_snp_baseline,
                                    (
                                        train_bfile, test_bfile,
                                        ot_ids, ov_ids,
                                        train_pheno, test_pheno,
                                        os.path.join(baseline_dir, "baseline_orig"),
                                        temp_snp_file
                                    ),
                                    {
                                        "_result_key": f"outer{outer_i}:baseline_orig",
                                        "description": "Original Density (GBLUP - Re-extracted)",
                                    },
                                    "Baseline Original Density (GBLUP - Re-extracted)"
                                )
                            except Exception as e:
                                LOGGER.log(f"Baseline-only: failed to launch baseline_orig on fold {outer_i}: {e}")
                        else:
                            LOGGER.log(f"Warning: No SNPs extracted from baseline bfile {baseline_bfile_arg}")
                    except Exception as e:
                        LOGGER.log(f"Baseline-only: baseline_orig failed on fold {outer_i}: {e}")

                if (not args.skip_hiblup_baseline) and args.baseline_full and fr.get("baseline_full_acc") is None:
                    try:
                        launch_baseline(
                            run_fixed_snp_baseline,
                            (
                                train_bfile, test_bfile,
                                ot_ids, ov_ids,
                                train_pheno, test_pheno,
                                os.path.join(baseline_dir, "baseline_full_input"),
                                None
                            ),
                            {
                                "_result_key": f"outer{outer_i}:baseline_full_input",
                                "description": "Full-Genome GBLUP (Input Data)",
                            },
                            "Baseline Full-Genome GBLUP (Input Data)"
                        )
                    except Exception as e:
                        LOGGER.log(f"Baseline-only: baseline_full_input failed on fold {outer_i}: {e}")

                    if args.baseline_train_bfile and args.baseline_test_bfile:
                        if os.path.exists(f"{args.baseline_train_bfile}.bed") and os.path.exists(f"{args.baseline_test_bfile}.bed"):
                            try:
                                launch_baseline(
                                    run_fixed_snp_baseline,
                                    (
                                        args.baseline_train_bfile, args.baseline_test_bfile,
                                        ot_ids, ov_ids,
                                        train_pheno, test_pheno,
                                        os.path.join(baseline_dir, "baseline_full_external"),
                                        None
                                    ),
                                    {
                                        "_result_key": f"outer{outer_i}:baseline_full_external",
                                        "description": f"Full-Genome GBLUP (External bfile: {os.path.basename(args.baseline_train_bfile)})",
                                    },
                                    f"Baseline Full-Genome GBLUP (External bfile: {os.path.basename(args.baseline_train_bfile)})"
                                )
                            except Exception as e:
                                LOGGER.log(f"Baseline-only: baseline_full_external failed on fold {outer_i}: {e}")

                if (not args.skip_hiblup_baseline) and (args.baseline_custom_snps or args.baseline_custom_bfile) and fr.get("baseline_custom_acc") is None:
                    final_snp_file = None
                    c_train = args.baseline_custom_bfile if args.baseline_custom_bfile else train_bfile
                    c_test = args.baseline_custom_bfile if args.baseline_custom_bfile else test_bfile
                    if args.baseline_custom_snps:
                        final_snp_file = normalize_snp_list_for_bfile(
                            args.baseline_custom_snps,
                            c_train,
                            os.path.join(baseline_dir, "baseline_custom")
                        )
                    try:
                        launch_baseline(
                            run_fixed_snp_baseline,
                            (
                                c_train, c_test,
                                ot_ids, ov_ids,
                                train_pheno, test_pheno,
                                os.path.join(baseline_dir, "baseline_custom"),
                                final_snp_file
                            ),
                            {
                                "_result_key": f"outer{outer_i}:baseline_custom",
                                "description": "Custom Low-Density Chip",
                            },
                            "Baseline Custom Low-Density Chip"
                        )
                    except Exception as e:
                        LOGGER.log(f"Baseline-only: baseline_custom failed on fold {outer_i}: {e}")

                if (not args.skip_hiblup_baseline) and args.baseline_custom_plus_selected and fr.get("baseline_custom_plus_selected_acc") is None:
                    try:
                        base_snps = []
                        if args.baseline_custom_snps:
                            base_file = normalize_snp_list_for_bfile(
                                args.baseline_custom_snps,
                                train_bfile,
                                os.path.join(baseline_dir, "baseline_custom_plus_selected")
                            )
                            base_snps = read_snp_list_file(base_file)
                        elif args.baseline_custom_bfile:
                            target_id_fmt = detect_id_format(train_bfile)
                            base_snps = get_snps_from_bim(args.baseline_custom_bfile, target_format=target_id_fmt)

                        selected_snp_file = os.path.join(run_dir, "folds", f"outer{outer_i}", "final", f"{run_id}.snps")
                        selected_snps = read_snp_list_file(selected_snp_file)
                        if not selected_snps:
                            final_dir = os.path.join(run_dir, "folds", f"outer{outer_i}", "final")
                            if os.path.isdir(final_dir):
                                cand = sorted(
                                    [os.path.join(final_dir, fn) for fn in os.listdir(final_dir) if fn.endswith(".snps")],
                                    key=lambda p: (run_id not in os.path.basename(p), len(os.path.basename(p))),
                                )
                                for fp in cand:
                                    selected_snps = read_snp_list_file(fp)
                                    if selected_snps:
                                        selected_snp_file = fp
                                        break

                        if not base_snps:
                            LOGGER.log(f"Baseline-only: baseline_custom_plus_selected skipped on fold {outer_i}: no custom SNPs")
                        elif not selected_snps:
                            LOGGER.log(f"Baseline-only: baseline_custom_plus_selected skipped on fold {outer_i}: missing selected SNPs file {selected_snp_file}")
                        else:
                            union = sorted(set(base_snps).union(selected_snps))
                            union_file = os.path.join(baseline_dir, f"baseline_custom_plus_selected_outer{outer_i}.snps")
                            pd.DataFrame(union, columns=["SNP"]).to_csv(union_file, index=False, header=False)
                            launch_baseline(
                                run_fixed_snp_baseline,
                                (
                                    train_bfile, test_bfile,
                                    ot_ids, ov_ids,
                                    train_pheno, test_pheno,
                                    os.path.join(baseline_dir, "baseline_custom_plus_selected"),
                                    union_file
                                ),
                                {
                                    "_result_key": f"outer{outer_i}:baseline_custom_plus_selected",
                                    "description": "Custom Chip SNPs ∪ Selected SNPs (Input bfile)",
                                },
                                "Baseline Custom+Selected (Input bfile)"
                            )
                    except Exception as e:
                        LOGGER.log(f"Baseline-only: baseline_custom_plus_selected failed on fold {outer_i}: {e}")

                if (not args.skip_hiblup_baseline) and args.baseline_commercial_snps:
                    c_train = train_bfile
                    c_test = test_bfile
                    chip_lists = args.baseline_commercial_snps.split(",")
                    for idx, chip_list in enumerate(chip_lists):
                        chip_list = chip_list.strip()
                        if not chip_list:
                            continue
                        comm_existing = fr.get("baseline_commercial") if isinstance(fr.get("baseline_commercial"), dict) else {}
                        if comm_existing.get(f"baseline_comm_{idx+1}") is not None:
                            continue
                        mapped_chip_list = normalize_snp_list_for_bfile(
                            chip_list,
                            c_train,
                            os.path.join(baseline_dir, f"baseline_comm_{idx+1}")
                        )
                        try:
                            launch_baseline(
                                run_fixed_snp_baseline,
                                (
                                    c_train, c_test,
                                    ot_ids, ov_ids,
                                    train_pheno, test_pheno,
                                    os.path.join(baseline_dir, f"baseline_comm_{idx+1}"),
                                    mapped_chip_list
                                ),
                                {
                                    "_result_key": f"outer{outer_i}:baseline_comm_{idx+1}",
                                    "description": f"Commercial Chip {idx+1}",
                                },
                                f"Baseline Commercial Chip {idx+1}"
                            )
                        except Exception as e:
                            LOGGER.log(f"Baseline-only: baseline_comm_{idx+1} failed on fold {outer_i}: {e}")

                target_id_fmt = detect_id_format(train_bfile)
                for idx, bl in enumerate(additional_baselines_list):
                    try:
                        bl_name = bl.get("name", f"Additional_Baseline_{idx}")
                        bl_bfile = bl.get("bfile")
                        bl_snps = bl.get("snp_list")
                        c_train = train_bfile
                        c_test = test_bfile
                        final_snp_file = None
                        if bl_snps and bl_snps not in ["all", "None", ""]:
                            final_snp_file = bl_snps
                        if not final_snp_file and bl_bfile and bl_bfile not in ["default", "None", ""]:
                            extracted_snps = get_snps_from_bim(bl_bfile, target_format=target_id_fmt)
                            if extracted_snps:
                                temp_snp_file = os.path.join(baseline_dir, f"baseline_add_{idx+1}.snps")
                                pd.DataFrame(extracted_snps, columns=["SNP"]).to_csv(temp_snp_file, index=False, header=False)
                                final_snp_file = temp_snp_file
                        if final_snp_file or bl_snps == "all":
                            launch_baseline(
                                run_fixed_snp_baseline,
                                (
                                    c_train, c_test,
                                    ot_ids, ov_ids,
                                    train_pheno, test_pheno,
                                    os.path.join(baseline_dir, f"baseline_add_{idx+1}_{bl_name.replace(' ', '_')}"),
                                    final_snp_file
                                ),
                                {
                                    "_result_key": f"outer{outer_i}:baseline_add_{idx+1}:{bl_name}",
                                    "description": f"Additional Baseline: {bl_name}",
                                },
                                f"Additional Baseline: {bl_name}"
                            )
                        else:
                            if LOGGER:
                                LOGGER.log(f"Skipping Additional Baseline {bl_name}: No SNPs found (bfile={bl_bfile}, snp_list={bl_snps})")
                    except Exception as e:
                        LOGGER.log(f"Baseline-only: additional baseline {bl.get('name')} failed on fold {outer_i}: {e}")

            try:
                for t in list(BASELINE_THREADS):
                    try:
                        t.join()
                    except Exception:
                        pass
            except Exception:
                pass

            try:
                with BASELINE_LOCK:
                    baseline_snapshot = dict(BASELINE_RESULTS)
            except Exception:
                baseline_snapshot = {}

            baseline_accuracies = []
            for fr in fold_records:
                try:
                    fold_no = int(fr.get("outer_fold"))
                except Exception:
                    continue
                m_full = baseline_snapshot.get(f"outer{fold_no}:baseline_full_input")
                if isinstance(m_full, dict):
                    fr["baseline_full_acc"] = safe_accuracy(m_full)
                    fr["baseline_full_r2"] = m_full.get("r2", None)
                    fr["baseline_full_spearman"] = m_full.get("spearman", None)
                if fr.get("baseline_full_acc") is not None:
                    try:
                        baseline_accuracies.append(float(fr["baseline_full_acc"]))
                    except Exception:
                        pass

                m_orig = baseline_snapshot.get(f"outer{fold_no}:baseline_orig")
                if isinstance(m_orig, dict):
                    fr["baseline_orig_acc"] = safe_accuracy(m_orig)

                m_custom = baseline_snapshot.get(f"outer{fold_no}:baseline_custom")
                if isinstance(m_custom, dict):
                    fr["baseline_custom_acc"] = safe_accuracy(m_custom)

                m_cps = baseline_snapshot.get(f"outer{fold_no}:baseline_custom_plus_selected")
                if isinstance(m_cps, dict):
                    fr["baseline_custom_plus_selected_acc"] = safe_accuracy(m_cps)

                comm = fr.get("baseline_commercial") if isinstance(fr.get("baseline_commercial"), dict) else {}
                for k, v in baseline_snapshot.items():
                    if not k.startswith(f"outer{fold_no}:baseline_comm_"):
                        continue
                    if isinstance(v, dict):
                        comm[k.split(":", 1)[1]] = safe_accuracy(v)
                if comm:
                    fr["baseline_commercial"] = comm
                    if fr.get("baseline_commercial_acc") is None:
                        try:
                            fr["baseline_commercial_acc"] = float(max(comm.values()))
                        except Exception:
                            pass

                adds = fr.get("baseline_additional") if isinstance(fr.get("baseline_additional"), dict) else {}
                for k, v in baseline_snapshot.items():
                    if not k.startswith(f"outer{fold_no}:baseline_add_"):
                        continue
                    if isinstance(v, dict):
                        adds[k.split(":", 1)[1]] = safe_accuracy(v)
                if adds:
                    fr["baseline_additional"] = adds

            try:
                os.makedirs(summary_dir, exist_ok=True)
            except Exception:
                pass
            trait = infer_trait_name_from_pheno(args.pheno or getattr(args, "train_pheno", None))
            with open(summary_txt, "w") as f:
                f.write(f"Trait: {trait}\n")
                f.write(f"RunID: {run_id}\n")
                f.write(f"OutDir: {run_dir}\n")
                f.write(f"Accuracies: {json.dumps(outer_accuracies)}\n")
                f.write(f"Mean: {mean_outer}\n")
                if baseline_accuracies:
                    mean_base = float(np.mean(baseline_accuracies))
                    f.write(f"Baseline_Accuracies: {json.dumps(baseline_accuracies)}\n")
                    f.write(f"Baseline_Mean: {mean_base}\n")
                f.write(f"Folds: {json.dumps(fold_records, ensure_ascii=False)}\n")
            try:
                with open(summary_json, "w") as f:
                    json.dump(
                        {
                            "trait": trait,
                            "run_id": run_id,
                            "out_dir": run_dir,
                            "outer_folds": int(folds_n),
                            "inner_folds": int(args.inner_cv_folds),
                            "accuracies": outer_accuracies,
                            "mean_accuracy": mean_outer,
                            "baseline_accuracies": baseline_accuracies,
                            "folds": fold_records,
                        },
                        f,
                        ensure_ascii=False,
                        indent=2,
                    )
            except Exception:
                pass
            return

        summary_dir = os.path.join(run_dir, "summary")
        summary_json = os.path.join(summary_dir, "cv_summary.json")
        fold_by_no = None
        if args.fold_index:
            fold_by_no = {}
            existing_summary = None
            if os.path.exists(summary_json):
                try:
                    with open(summary_json, "r") as f:
                        existing_summary = json.load(f)
                except Exception:
                    existing_summary = None
            if isinstance(existing_summary, dict):
                try:
                    existing_folds = list(existing_summary.get("folds") or [])
                except Exception:
                    existing_folds = []
                for fr in existing_folds:
                    if not isinstance(fr, dict):
                        continue
                    try:
                        fold_no = int(fr.get("outer_fold"))
                    except Exception:
                        continue
                    fold_by_no[fold_no] = fr

        outer_kf = KFold(n_splits=args.folds, shuffle=True, random_state=42)
        X_indices = np.arange(len(valid_train))
        outer_accuracies = []
        baseline_accuracies = []
        fold_records = []
        
        for outer_i, (tr_idx, te_idx) in enumerate(outer_kf.split(X_indices), start=1):
            if args.fold_index and outer_i != args.fold_index:
                continue

            outer_train = valid_train.iloc[tr_idx]
            outer_test = valid_train.iloc[te_idx]
            outer_dir = os.path.join(run_dir, "folds", f"outer{outer_i}")
            try:
                os.makedirs(outer_dir, exist_ok=True)
                os.makedirs(os.path.join(outer_dir, "ids"), exist_ok=True)
            except Exception:
                pass
            fold_record = None
            if fold_by_no is not None:
                try:
                    existing_fr = fold_by_no.get(int(outer_i))
                except Exception:
                    existing_fr = None
                if isinstance(existing_fr, dict):
                    fold_record = dict(existing_fr)
            if not isinstance(fold_record, dict):
                fold_record = {"outer_fold": int(outer_i)}
            fold_record["train_n"] = int(len(outer_train))
            fold_record["test_n"] = int(len(outer_test))
            
            LOGGER.log(f"=========================================================")
            LOGGER.log(f"                OUTER FOLD {outer_i}/{args.folds}")
            LOGGER.log(f"=========================================================")
            LOGGER.log(f"Training Set Size: {len(outer_train)} individuals")
            LOGGER.log(f"Testing Set Size:  {len(outer_test)} individuals")
            
            # Write outer ID files
            ot_ids = os.path.join(outer_dir, "ids", "train.ids")
            ov_ids = os.path.join(outer_dir, "ids", "test.ids")
            outer_train[["FID","IID"]].to_csv(ot_ids, sep="\t", index=False, header=False)
            outer_test[["FID","IID"]].to_csv(ov_ids, sep="\t", index=False, header=False)

            if args.fast_mode:
                p_space = [1e-5, 1e-4, 1e-3, 0.01, 0.05]
                r2_space = [0.1, 0.2, 0.5]
                kb_space = [250, 500]
            else:
                p_space = [1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
                r2_space = [0.1, 0.2, 0.5, 0.8]
                kb_space = [250, 500, 1000]
            max_snps_fold = int(args.max_snps) if args.max_snps else None

            selection_method_fold = args.selection_method
            need_background = bool(args.auto_search_space) or (getattr(args, "architecture_preset", "auto") != "auto") or (args.selection_method == "auto") or bool(prior)
            if need_background:
                prior_preset = None
                if prior.get("preset") in {"polygenic", "oligogenic", "hybrid"}:
                    prior_preset = prior["preset"]
                elif prior.get("architecture_preset") in {"polygenic", "oligogenic", "hybrid"}:
                    prior_preset = prior["architecture_preset"]

                preset = args.architecture_preset if args.architecture_preset != "auto" else prior_preset
                diag = None
                if not preset:
                    diag_file = None
                    diag_is_fastgwa = False
                    if external_qassoc and os.path.exists(external_qassoc):
                        diag_file = external_qassoc
                        diag_is_fastgwa = False
                    else:
                        try:
                            os.makedirs(os.path.join(outer_dir, "diag"), exist_ok=True)
                        except Exception:
                            pass
                        diag_prefix = os.path.join(outer_dir, "diag", run_id)
                        diag_file, diag_is_fastgwa, _ = run_gwas_only(
                            train_bfile, train_pheno, ot_ids, diag_prefix,
                            args.plink_threads, use_fastgwa=True, background=False, reuse=args.reuse_gwas
                        )
                    if diag_file and os.path.exists(diag_file):
                        diag = gwas_diagnostics(diag_file, is_fastgwa=diag_is_fastgwa)
                        preset = infer_architecture_preset(diag)
                    else:
                        preset = "hybrid"

                space = biased_search_space(
                    preset=preset,
                    fast_mode=bool(args.fast_mode),
                    base_max_snps=int(args.max_snps),
                    base_clump_cap=int(args.clump_max_summary_snps),
                )
                space = apply_prior_to_space(space, prior)
                p_space = space["p_space"]
                r2_space = space["r2_space"]
                kb_space = space["kb_space"]
                max_snps_fold = space["max_snps"]
                try:
                    CLUMP_MAX_SUMMARY_SNPS = int(space["clump_cap"])
                except Exception:
                    pass

                if args.selection_method == "auto":
                    selection_method_fold = choose_selection_method(space["preset"], prior=prior)

                if LOGGER:
                    if diag and diag.get("n"):
                        LOGGER.log(
                            f"Architecture preset: {space['preset']} (N={diag['n']}, p<1e-4={diag['counts'].get(1e-4,0)}, lambdaGC={diag.get('lambda_gc')}); selection_method={selection_method_fold}"
                        )
                    else:
                        LOGGER.log(f"Architecture preset: {space['preset']}; selection_method={selection_method_fold}")
            elif args.selection_method == "auto":
                selection_method_fold = "hybrid"
            selection_method_candidates_fold = None
            if args.optimize_selection_method and args.selection_method == "auto":
                preset_for_mode = None
                try:
                    preset_for_mode = (space.get("preset") if isinstance(space, dict) else None) if need_background else None
                except Exception:
                    preset_for_mode = None
                if not preset_for_mode:
                    try:
                        preset_for_mode = (args.architecture_preset if getattr(args, "architecture_preset", "auto") != "auto" else None) or "hybrid"
                    except Exception:
                        preset_for_mode = "hybrid"
                preset_for_mode = str(preset_for_mode).lower()
                if preset_for_mode == "oligogenic":
                    selection_method_candidates_fold = ["clump", "hybrid"]
                elif preset_for_mode == "polygenic":
                    selection_method_candidates_fold = ["prune", "hybrid"]
                else:
                    selection_method_candidates_fold = ["clump", "prune", "hybrid"]
            
            # ---------------------------------------------------------
            # Pre-calculate Inner CV Splits and Launch GWAS (Parallel)
            # ---------------------------------------------------------
            inner_cv_folds = args.inner_cv_folds
            inner_kf = KFold(n_splits=inner_cv_folds, shuffle=True, random_state=42)
            inner_gwas_processes = []
            
            LOGGER.log(f"Pre-calculating {inner_cv_folds} Inner CV folds and launching GWAS in background...")
            
            # Reconstruct indices for inner split relative to outer_train
            inner_split_indices = list(inner_kf.split(np.arange(len(outer_train))))
            
            for i, (it_idx, iv_idx) in enumerate(inner_split_indices):
                # Prepare IDs
                inner_train = outer_train.iloc[it_idx]
                inner_val = outer_train.iloc[iv_idx]
                inner_dir = os.path.join(outer_dir, "inner", f"fold{i+1}")
                try:
                    os.makedirs(os.path.join(inner_dir, "ids"), exist_ok=True)
                    os.makedirs(os.path.join(inner_dir, "gwas"), exist_ok=True)
                except Exception:
                    pass
                it_ids = os.path.join(inner_dir, "ids", "train.ids")
                iv_ids = os.path.join(inner_dir, "ids", "val.ids")
                
                # Always write (or check exist)
                if not os.path.exists(it_ids):
                     inner_train[["FID","IID"]].to_csv(it_ids, sep="\t", index=False, header=False)
                if not os.path.exists(iv_ids):
                     inner_val[["FID","IID"]].to_csv(iv_ids, sep="\t", index=False, header=False)
                     
                # Prepare/Reuse GWAS if not external
                if not external_qassoc:
                    gwas_shared_prefix = os.path.join(inner_dir, "gwas", run_id)
                    use_fastgwa = True # Enable fastGWA by default for speed
                    
                    target_gwas_file, is_fastgwa_out, p = run_gwas_only(
                        train_bfile, train_pheno, it_ids, gwas_shared_prefix, 
                        args.plink_threads, use_fastgwa=use_fastgwa, background=True, reuse=args.reuse_gwas
                    )
                    
                    if p:
                        inner_gwas_processes.append((i, p, target_gwas_file, is_fastgwa_out))
                        LOGGER.log(f"Launching Background {'PLINK2' if is_fastgwa_out else 'PLINK1.9'} GWAS for Inner Fold {i+1}...")
                    elif os.path.exists(target_gwas_file):
                         LOGGER.log(f"Reusing Existing GWAS for Inner Fold {i+1} ({'PLINK2' if is_fastgwa_out else 'PLINK1.9'}): {target_gwas_file}")
                    else:
                         LOGGER.log(f"Failed to launch or find GWAS for Inner Fold {i+1}")

            # Prepare external GWAS qassoc if provided (reused)
            # (Logic handled above, external_qassoc is already set)
            
            # ---------------------------------------------------------
            # Baselines
            # ---------------------------------------------------------
            
            # 1. Existing "Original Density" Baseline (if provided via --baseline_bfile)
            baseline_bfile_arg = args.baseline_train_bfile if args.baseline_train_bfile else args.baseline_bfile
            baseline_dir = os.path.join(outer_dir, "baselines")
            try:
                os.makedirs(baseline_dir, exist_ok=True)
            except Exception:
                pass
            
            if not args.skip_hiblup_baseline and baseline_bfile_arg:
                try:
                    # Logic: Use study population (train_bfile), extract SNPs from baseline_bfile_arg
                    target_id_fmt = detect_id_format(train_bfile)
                    extracted_snps = get_snps_from_bim(baseline_bfile_arg, target_format=target_id_fmt)
                    
                    if extracted_snps:
                        temp_snp_file = os.path.join(baseline_dir, "baseline_orig.snps")
                        pd.DataFrame(extracted_snps, columns=["SNP"]).to_csv(temp_snp_file, index=False, header=False)
                        try:
                            launch_baseline(
                                run_fixed_snp_baseline,
                                (
                                    train_bfile, test_bfile,
                                    ot_ids, ov_ids,
                                    train_pheno, test_pheno,
                                    os.path.join(baseline_dir, "baseline_orig"),
                                    temp_snp_file
                                ),
                                {
                                    "_result_key": f"outer{outer_i}:baseline_orig",
                                    "description": "Original Density (GBLUP - Re-extracted)",
                                },
                                "Baseline Original Density (GBLUP - Re-extracted)"
                            )
                        except Exception as e:
                            LOGGER.log(f"Failed to launch baseline_orig: {e}")
                    else:
                        LOGGER.log(f"Warning: No SNPs extracted from baseline bfile {baseline_bfile_arg}")
                        
                except Exception as e:
                    LOGGER.log(f"Baseline GBLUP failed on fold {outer_i}: {e}")

            # 2. Full-Genome GBLUP (on Input Data)
            if (not args.skip_hiblup_baseline) and args.baseline_full:
                try:
                    launch_baseline(
                        run_fixed_snp_baseline,
                        (
                            train_bfile, test_bfile,
                            ot_ids, ov_ids,
                            train_pheno, test_pheno,
                            os.path.join(baseline_dir, "baseline_full_input"),
                            None
                        ),
                        {
                            "_result_key": f"outer{outer_i}:baseline_full_input",
                            "description": "Full-Genome GBLUP (Input Data)",
                        },
                        "Baseline Full-Genome GBLUP (Input Data)"
                    )
                except Exception as e:
                    LOGGER.log(f"Full-Genome Baseline failed on fold {outer_i}: {e}")

                if args.baseline_train_bfile and args.baseline_test_bfile:
                    if os.path.exists(f"{args.baseline_train_bfile}.bed") and os.path.exists(f"{args.baseline_test_bfile}.bed"):
                        try:
                            launch_baseline(
                                run_fixed_snp_baseline,
                                (
                                    args.baseline_train_bfile, args.baseline_test_bfile,
                                    ot_ids, ov_ids,
                                    train_pheno, test_pheno,
                                    os.path.join(baseline_dir, "baseline_full_external"),
                                    None
                                ),
                                {
                                    "_result_key": f"outer{outer_i}:baseline_full_external",
                                    "description": f"Full-Genome GBLUP (External bfile: {os.path.basename(args.baseline_train_bfile)})",
                                },
                                f"Baseline Full-Genome GBLUP (External bfile: {os.path.basename(args.baseline_train_bfile)})"
                            )
                        except Exception as e:
                            LOGGER.log(f"Full-Genome Baseline (external bfile) failed on fold {outer_i}: {e}")

            # 3. Custom Low-Density Chip (No Imputation)
            if (not args.skip_hiblup_baseline) and (args.baseline_custom_snps or args.baseline_custom_bfile):
                final_snp_file = None

                c_train = args.baseline_custom_bfile if args.baseline_custom_bfile else train_bfile
                c_test = args.baseline_custom_bfile if args.baseline_custom_bfile else test_bfile

                if args.baseline_custom_snps:
                    final_snp_file = normalize_snp_list_for_bfile(
                        args.baseline_custom_snps,
                        c_train,
                        os.path.join(baseline_dir, "baseline_custom")
                    )
                
                try:
                    launch_baseline(
                        run_fixed_snp_baseline,
                        (
                            c_train, c_test,
                            ot_ids, ov_ids,
                            train_pheno, test_pheno,
                            os.path.join(baseline_dir, "baseline_custom"),
                            final_snp_file
                        ),
                        {
                            "_result_key": f"outer{outer_i}:baseline_custom",
                            "description": "Custom Low-Density Chip",
                        },
                        "Baseline Custom Low-Density Chip"
                    )
                except Exception as e:
                    LOGGER.log(f"Custom Chip Baseline failed on fold {outer_i}: {e}")

            # 4. Commercial Chips
            if (not args.skip_hiblup_baseline) and args.baseline_commercial_snps:
                # FORCE use of input data
                c_train = train_bfile
                c_test = test_bfile
                
                # Handle comma-separated list
                chip_lists = args.baseline_commercial_snps.split(",")
                for idx, chip_list in enumerate(chip_lists):
                    chip_list = chip_list.strip()
                    if not chip_list: continue
                    mapped_chip_list = normalize_snp_list_for_bfile(
                        chip_list,
                        c_train,
                        os.path.join(baseline_dir, f"baseline_comm_{idx+1}")
                    )
                    
                    try:
                        launch_baseline(
                            run_fixed_snp_baseline,
                            (
                                c_train, c_test,
                                ot_ids, ov_ids,
                                train_pheno, test_pheno,
                                os.path.join(baseline_dir, f"baseline_comm_{idx+1}"),
                                mapped_chip_list
                            ),
                            {
                                "_result_key": f"outer{outer_i}:baseline_comm_{idx+1}",
                                "description": f"Commercial Chip {idx+1}",
                            },
                            f"Baseline Commercial Chip {idx+1}"
                        )
                    except Exception as e:
                        LOGGER.log(f"Commercial Chip Baseline {idx+1} failed on fold {outer_i}: {e}")
            
            # 5. Additional Baselines
            
            # Detect Target Format
            target_id_fmt = detect_id_format(train_bfile)
            
            for idx, bl in enumerate(additional_baselines_list):
                try:
                    bl_name = bl.get("name", f"Additional_Baseline_{idx}")
                    bl_bfile = bl.get("bfile")
                    bl_snps = bl.get("snp_list")
                    
                    # Force data source to be the current population
                    c_train = train_bfile
                    c_test = test_bfile
                    
                    final_snp_file = None
                    
                    # Handle SNP source
                    # If snp_list is provided, use it
                    if bl_snps and bl_snps not in ["all", "None", ""]:
                        final_snp_file = bl_snps
                    
                    # If bfile is provided, extract SNPs from it if snp_list is not provided
                    if not final_snp_file and bl_bfile and bl_bfile not in ["default", "None", ""]:
                        extracted_snps = get_snps_from_bim(bl_bfile, target_format=target_id_fmt)
                        if extracted_snps:
                            temp_snp_file = os.path.join(baseline_dir, f"baseline_add_{idx+1}.snps")
                            pd.DataFrame(extracted_snps, columns=["SNP"]).to_csv(temp_snp_file, index=False, header=False)
                            final_snp_file = temp_snp_file
                    
                    if final_snp_file or bl_snps == "all":
                         launch_baseline(
                            run_fixed_snp_baseline,
                            (
                                c_train, c_test,
                                ot_ids, ov_ids,
                                train_pheno, test_pheno,
                                os.path.join(baseline_dir, f"baseline_add_{idx+1}_{bl_name.replace(' ', '_')}"),
                                final_snp_file
                            ),
                            {
                                "_result_key": f"outer{outer_i}:baseline_add_{idx+1}:{bl_name}",
                                "description": f"Additional Baseline: {bl_name}",
                            },
                            f"Additional Baseline: {bl_name}"
                         )
                    else:
                        if LOGGER: LOGGER.log(f"Skipping Additional Baseline {bl_name}: No SNPs found (bfile={bl_bfile}, snp_list={bl_snps})")

                except Exception as e:
                    LOGGER.log(f"Additional Baseline {bl.get('name')} failed on fold {outer_i}: {e}")

            # ---------------------------------------------------------
            # Wait for Background GWAS
            # ---------------------------------------------------------
            if inner_gwas_processes:
                LOGGER.log(f"Waiting for {len(inner_gwas_processes)} background GWAS jobs to complete...")
                for i, p, fpath, is_plink2_out in inner_gwas_processes:
                    p.wait()
                    if p.returncode != 0:
                        LOGGER.log(f"Error: Background GWAS failed for Fold {i+1}")
                    else:
                        LOGGER.log(f"Background GWAS completed for Fold {i+1}")
                        if is_plink2_out:
                            files = glob.glob(fpath)
                            if files:
                                actual_file = files[0]
                                std_name = fpath.replace("*.glm.linear", "fastGWA")
                                if actual_file != std_name:
                                    try:
                                        if os.path.exists(std_name):
                                            os.remove(std_name)
                                        os.rename(actual_file, std_name)
                                        LOGGER.log(f"Renamed {actual_file} to {std_name}")
                                    except Exception as e:
                                        LOGGER.log(f"Error renaming PLINK2 output: {e}")
                                conv_out = f"{std_name}.plink_fmt.temp"
                                need_conv = True
                                if os.path.exists(conv_out):
                                    try:
                                        with open(conv_out, "r") as f:
                                            _ = f.readline()
                                            _ = f.readline()
                                        need_conv = True if (_ is None or _ == "") else False
                                    except Exception:
                                        need_conv = True
                                if need_conv:
                                    cmd_conv = (
                                        "awk 'BEGIN{OFS=\"\\t\"} "
                                        "NR==1{for(i=1;i<=NF;i++){if($i==\"SNP\"||$i==\"ID\")s=i; if($i==\"P\")p=i; if($i==\"TEST\")t=i;} "
                                        "print \"SNP\",\"P\"; next} "
                                        "{if(s>0 && p>0){if(t==0 || $t==\"ADD\"){if($p!=\"NA\" && $p!=\"\"){print $s,$p}}}}' "
                                        f"{std_name} > {conv_out}"
                                    )
                                    if run_command(cmd_conv):
                                        LOGGER.log(f"Created lightweight GWAS file {conv_out}")
                                    else:
                                        LOGGER.log(f"Failed to create lightweight GWAS file {conv_out}")
                            else:
                                LOGGER.log(f"Error: PLINK2 output not found for pattern {fpath}")
                inner_gwas_processes = []

            if selection_method_fold == "clump_cached":
                try:
                    cache_pmax = float(max(p_space))
                    cache_kb = int(max(kb_space))
                    cache_min_r2 = float(min(r2_space))
                    LOGGER.log(f"Prewarming cached clump LD (pmax={cache_pmax}, kb={cache_kb}, min_r2={cache_min_r2}, max_snps={CLUMP_CACHED_MAX_SNPS})")
                    for i in range(args.inner_cv_folds):
                        inner_dir = os.path.join(outer_dir, "inner", f"fold{i+1}")
                        gwas_shared_prefix = os.path.join(inner_dir, "gwas", run_id)
                        gwas_assoc_file_plink = f"{gwas_shared_prefix}.assoc.linear"
                        gwas_assoc_file_fastgwa = f"{gwas_shared_prefix}.fastGWA"
                        gwas_assoc_file_fastgwa_conv = f"{gwas_shared_prefix}.fastGWA.plink_fmt.temp"
                        gwas_clump_file = None
                        is_fastgwa_out = False
                        if os.path.exists(gwas_assoc_file_fastgwa_conv):
                            gwas_clump_file = gwas_assoc_file_fastgwa_conv
                            is_fastgwa_out = False
                        elif os.path.exists(gwas_assoc_file_fastgwa):
                            gwas_clump_file = gwas_assoc_file_fastgwa
                            is_fastgwa_out = True
                        elif os.path.exists(gwas_assoc_file_plink):
                            gwas_clump_file = gwas_assoc_file_plink
                            is_fastgwa_out = False
                        else:
                            plink2_files = glob.glob(f"{gwas_shared_prefix}.*.glm.linear")
                            if plink2_files:
                                gwas_clump_file = plink2_files[0]
                                is_fastgwa_out = True
                        if not gwas_clump_file or not os.path.exists(gwas_clump_file):
                            continue
                        _ = _get_or_build_clump_ld_cache(
                            bfile=train_bfile,
                            gwas_file=gwas_clump_file,
                            is_fastgwa=is_fastgwa_out,
                            out_prefix=gwas_shared_prefix,
                            p_max=cache_pmax,
                            max_kb=cache_kb,
                            min_r2=cache_min_r2,
                            plink_threads=PLINK_THREADS,
                        )
                except Exception as e:
                    LOGGER.log(f"Cached clump LD prewarm failed: {e}")

            # Inner CV Objective on outer_train
            def objective(trial):
                # ... (rest of objective function)
                # Define method search space if method="auto"
                if args.method == "auto":
                    if args.proxy_metric:
                        method = "rrblup"
                    else:
                        method = trial.suggest_categorical("method", ["rrblup", "adaptive", "lightgbm"])
                else:
                    method = args.method

                p_thresh = trial.suggest_categorical("p_thresh", p_space)
                clump_r2 = trial.suggest_categorical("clump_r2", r2_space)
                clump_kb = trial.suggest_categorical("clump_kb", kb_space)
                info_thresh = trial.suggest_categorical("info_thresh", [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]) if args.info_file else 0.0
                
                n_folds_cv = args.inner_cv_folds
                kf = KFold(n_splits=n_folds_cv, shuffle=True, random_state=42)
                accuracies = []
                
                idx = np.arange(len(outer_train))
                for i, (it_idx, iv_idx) in enumerate(kf.split(idx)):
                    LOGGER.log(f"[Trial {trial.number}] --- Internal CV Fold {i+1}/{n_folds_cv} ---")
                    inner_train = outer_train.iloc[it_idx]
                    inner_val = outer_train.iloc[iv_idx]
                    
                    # Use shared IDs for Inner CV to avoid duplication across trials
                    inner_dir = os.path.join(outer_dir, "inner", f"fold{i+1}")
                    try:
                        os.makedirs(os.path.join(inner_dir, "ids"), exist_ok=True)
                        os.makedirs(os.path.join(inner_dir, "gwas"), exist_ok=True)
                    except Exception:
                        pass
                    it_ids = os.path.join(inner_dir, "ids", "train.ids")
                    iv_ids = os.path.join(inner_dir, "ids", "val.ids")
                    
                    if not os.path.exists(it_ids):
                        inner_train[["FID","IID"]].to_csv(it_ids, sep="\t", index=False, header=False)
                    if not os.path.exists(iv_ids):
                        inner_val[["FID","IID"]].to_csv(iv_ids, sep="\t", index=False, header=False)
                    
                    # Prepare Extract File (Info Score)
                    extract_file = ensure_info_extract(args.info_file, info_thresh, args.out, target_bfile=train_bfile) if (args.info_file and info_thresh > 0) else None
                    
                    # Trial-specific prefix for Clumping output
                    trial_dir = os.path.join(outer_dir, "trials", f"trial{trial.number}", f"fold{i+1}")
                    try:
                        os.makedirs(trial_dir, exist_ok=True)
                    except Exception:
                        pass
                    gwas_prefix = os.path.join(trial_dir, run_id)
                    
                    # Shared GWAS prefix for reuse
                    gwas_shared_prefix = os.path.join(inner_dir, "gwas", run_id)
                    gwas_assoc_file_plink = f"{gwas_shared_prefix}.assoc.linear"
                    gwas_assoc_file_fastgwa = f"{gwas_shared_prefix}.fastGWA"
                    gwas_assoc_file_fastgwa_conv = f"{gwas_shared_prefix}.fastGWA.plink_fmt.temp"
                    
                    # Decide which file to use for SNP selection (clump/prune)
                    qassoc_source = None
                    if os.path.exists(gwas_assoc_file_fastgwa_conv):
                        gwas_clump_file = gwas_assoc_file_fastgwa_conv
                        is_fastgwa_out = False
                        qassoc_source = gwas_assoc_file_fastgwa if os.path.exists(gwas_assoc_file_fastgwa) else gwas_clump_file
                    elif os.path.exists(gwas_assoc_file_fastgwa):
                        pfile = ensure_two_col_gwas_file(gwas_assoc_file_fastgwa, is_fastgwa=True)
                        gwas_clump_file = pfile if pfile else gwas_assoc_file_fastgwa
                        is_fastgwa_out = False if pfile else True
                        qassoc_source = gwas_assoc_file_fastgwa
                    elif os.path.exists(gwas_assoc_file_plink):
                        gwas_clump_file = gwas_assoc_file_plink
                        is_fastgwa_out = False
                        qassoc_source = gwas_assoc_file_plink
                    else:
                        plink2_files = glob.glob(f"{gwas_shared_prefix}.*.glm.linear")
                        if plink2_files:
                            qassoc_source = plink2_files[0]
                            pfile = ensure_two_col_gwas_file(qassoc_source, is_fastgwa=True)
                            gwas_clump_file = pfile if pfile else qassoc_source
                            is_fastgwa_out = False if pfile else True
                        else:
                            LOGGER.log(f"Error: Shared GWAS file missing for fold {i+1} ({gwas_assoc_file_plink} or {gwas_assoc_file_fastgwa} or {gwas_shared_prefix}.*.glm.linear)")
                            accuracies.append(0.0)
                            continue
                        
                    LOGGER.log(f"Reusing Shared GWAS for Inner Fold {i+1} ({'PLINK2' if is_fastgwa_out else 'PLINK1.9'})...")

                    selection_method_trial = selection_method_fold
                    if selection_method_candidates_fold:
                        selection_method_trial = trial.suggest_categorical("selection_method", selection_method_candidates_fold)
                    snps = select_snps_from_gwas(
                        bfile=train_bfile,
                        gwas_file=gwas_clump_file,
                        keep_ids_file=it_ids,
                        out_prefix=gwas_prefix,
                        p_thresh=p_thresh,
                        ld_r2=clump_r2,
                        window_kb=clump_kb,
                        plink_threads=PLINK_THREADS,
                        extract_file=extract_file,
                        is_fastgwa=is_fastgwa_out,
                        selection_method=selection_method_trial,
                        max_snps=max_snps_fold,
                    )
                    LOGGER.log(f"Selected {len(snps)} SNPs.")

                    extract_file = extract_file
                    
                    if len(snps) == 0:
                        LOGGER.log("No SNPs selected, skipping training.")
                        accuracies.append(0.0)
                        running_mean = float(np.mean(accuracies)) if accuracies else 0.0
                        trial.report(running_mean, step=i)
                        if trial.should_prune():
                            LOGGER.log(f"Pruned Trial {trial.number} at inner fold {i+1} (mean_acc={running_mean:.4f}).")
                            raise optuna.TrialPruned()
                        continue
                    
                    snp_file = f"{gwas_prefix}.snps"
                    pd.DataFrame(snps, columns=["SNP"]).to_csv(snp_file, index=False, header=False)
                    
                    if args.proxy_metric:
                        if external_qassoc:
                            qassoc_for_proxy = external_qassoc
                        else:
                            # Both SBLUP and PRS need full summary columns (A1/BETA)
                            # The .temp file might only have SNP/P for clumping.
                            qassoc_for_proxy = qassoc_source if qassoc_source else gwas_clump_file
                        
                        if args.proxy_method == "sblup":
                            metrics = calculate_sblup_proxy(
                                train_bfile, train_bfile,
                                it_ids, iv_ids,
                                train_pheno,
                                snp_file,
                                qassoc_for_proxy,
                                f"{gwas_prefix}.proxy",
                                gcta_path=args.gcta_path
                            )
                            acc = metrics["accuracy"]
                            LOGGER.log(f"Fold {i+1} Proxy Accuracy (SBLUP): {acc:.4f}")
                        else:
                            metrics = calculate_prs_proxy(
                                train_bfile, train_bfile,
                                it_ids, iv_ids,
                                train_pheno,
                                snp_file,
                                qassoc_for_proxy,
                                f"{gwas_prefix}.proxy"
                            )
                            acc = metrics["accuracy"]
                            LOGGER.log(f"Fold {i+1} Proxy Accuracy (PRS): {acc:.4f}")
                    else:
                        LOGGER.log("Training model...")
                        metrics = evaluate_performance(
                            train_bfile, train_bfile, train_pheno, train_pheno,
                            it_ids, iv_ids, snp_file, f"{gwas_prefix}.res",
                            method=method, libpath=args.libpath
                        )
                        acc_raw = metrics.get("accuracy", float("nan"))
                        acc = safe_accuracy(metrics)
                        if isinstance(acc_raw, float) and np.isnan(acc_raw):
                            LOGGER.log(f"Fold {i+1} Accuracy: NA (Method: {method})")
                        else:
                            LOGGER.log(f"Fold {i+1} Accuracy: {acc:.4f} (Method: {method})")
                    
                    accuracies.append(acc)
                    running_mean = float(np.mean(accuracies)) if accuracies else 0.0
                    trial.report(running_mean, step=i)
                    if trial.should_prune():
                        LOGGER.log(f"Pruned Trial {trial.number} at inner fold {i+1} (mean_acc={running_mean:.4f}).")
                        raise optuna.TrialPruned()
                
                mean_acc = float(np.mean(accuracies)) if len(accuracies) else 0.0
                LOGGER.log(f"Trial {trial.number} Mean Accuracy: {mean_acc:.4f}")
                return mean_acc
            
            pruner = optuna.pruners.MedianPruner(
                n_startup_trials=max(0, int(args.optuna_startup_trials)),
                n_warmup_steps=max(1, int(args.inner_cv_folds) // 2),
            ) if args.proxy_metric else None
            callbacks = [make_optuna_early_stop_callback(args.optuna_patience, args.optuna_min_delta)]
            study = optuna.create_study(
                direction="maximize",
                sampler=optuna.samplers.TPESampler(seed=args.optuna_seed),
                pruner=pruner,
            )
            study.optimize(objective, n_trials=args.trials, n_jobs=args.n_jobs, callbacks=callbacks)
            best_params = study.best_params
            best_params_proxy = dict(best_params) if best_params else {}

            if args.proxy_metric and args.two_stage_opt:
                try:
                    try:
                        topk = int(getattr(args, "two_stage_topk", 8))
                    except Exception:
                        topk = 8
                    topk = max(1, topk)
                    completed = [
                        t for t in study.trials
                        if getattr(t, "value", None) is not None and getattr(t, "state", None) == optuna.trial.TrialState.COMPLETE
                    ]
                    completed.sort(key=lambda t: float(t.value), reverse=True)
                    candidates = completed[:topk]
                    if candidates:
                        LOGGER.log(f"Two-stage: reevaluating top {len(candidates)} proxy trials with full metric...")

                    def _full_eval_params(params: dict, trial_no: int, rank: int) -> float:
                        p_thresh = params.get("p_thresh")
                        clump_r2 = params.get("clump_r2")
                        clump_kb = params.get("clump_kb")
                        try:
                            info_thresh = float(params.get("info_thresh", 0.0))
                        except Exception:
                            info_thresh = 0.0
                        method_full = resolve_method(args.method, params, proxy_metric=False)
                        n_folds_cv = args.inner_cv_folds
                        kf = KFold(n_splits=n_folds_cv, shuffle=True, random_state=42)
                        accuracies = []
                        idx = np.arange(len(outer_train))

                        for i, (it_idx, iv_idx) in enumerate(kf.split(idx)):
                            inner_train = outer_train.iloc[it_idx]
                            inner_val = outer_train.iloc[iv_idx]
                            inner_dir = os.path.join(outer_dir, "inner", f"fold{i+1}")
                            it_ids = os.path.join(inner_dir, "ids", "train.ids")
                            iv_ids = os.path.join(inner_dir, "ids", "val.ids")
                            if not os.path.exists(it_ids):
                                inner_train[["FID", "IID"]].to_csv(it_ids, sep="\t", index=False, header=False)
                            if not os.path.exists(iv_ids):
                                inner_val[["FID", "IID"]].to_csv(iv_ids, sep="\t", index=False, header=False)

                            extract_file = ensure_info_extract(args.info_file, info_thresh, args.out, target_bfile=train_bfile) if (args.info_file and info_thresh > 0) else None
                            two_stage_dir = os.path.join(outer_dir, "two_stage", f"trial{trial_no}", f"rank{rank}", f"fold{i+1}")
                            try:
                                os.makedirs(two_stage_dir, exist_ok=True)
                            except Exception:
                                pass
                            gwas_prefix = os.path.join(two_stage_dir, run_id)

                            gwas_shared_prefix = os.path.join(inner_dir, "gwas", run_id)
                            gwas_assoc_file_plink = f"{gwas_shared_prefix}.assoc.linear"
                            gwas_assoc_file_fastgwa = f"{gwas_shared_prefix}.fastGWA"
                            gwas_assoc_file_fastgwa_conv = f"{gwas_shared_prefix}.fastGWA.plink_fmt.temp"

                            if os.path.exists(gwas_assoc_file_fastgwa_conv):
                                gwas_clump_file = gwas_assoc_file_fastgwa_conv
                                is_fastgwa_out = False
                            elif os.path.exists(gwas_assoc_file_fastgwa):
                                pfile = ensure_two_col_gwas_file(gwas_assoc_file_fastgwa, is_fastgwa=True)
                                gwas_clump_file = pfile if pfile else gwas_assoc_file_fastgwa
                                is_fastgwa_out = False if pfile else True
                            elif os.path.exists(gwas_assoc_file_plink):
                                gwas_clump_file = gwas_assoc_file_plink
                                is_fastgwa_out = False
                            else:
                                plink2_files = glob.glob(f"{gwas_shared_prefix}.*.glm.linear")
                                if plink2_files:
                                    qassoc_source = plink2_files[0]
                                    pfile = ensure_two_col_gwas_file(qassoc_source, is_fastgwa=True)
                                    gwas_clump_file = pfile if pfile else qassoc_source
                                    is_fastgwa_out = False if pfile else True
                                else:
                                    accuracies.append(0.0)
                                    continue

                            snps = select_snps_from_gwas(
                                bfile=train_bfile,
                                gwas_file=gwas_clump_file,
                                keep_ids_file=it_ids,
                                out_prefix=gwas_prefix,
                                p_thresh=p_thresh,
                                ld_r2=clump_r2,
                                window_kb=clump_kb,
                                plink_threads=PLINK_THREADS,
                                extract_file=extract_file,
                                is_fastgwa=is_fastgwa_out,
                                selection_method=(
                                    str(params.get("selection_method"))
                                    if isinstance(params.get("selection_method"), str) and params.get("selection_method") in {"clump", "prune", "hybrid", "clump_cached"}
                                    else selection_method_fold
                                ),
                                max_snps=max_snps_fold,
                            )
                            if len(snps) == 0:
                                accuracies.append(0.0)
                                continue

                            snp_file = f"{gwas_prefix}.snps"
                            pd.DataFrame(snps, columns=["SNP"]).to_csv(snp_file, index=False, header=False)
                            metrics = evaluate_performance(
                                train_bfile, train_bfile, train_pheno, train_pheno,
                                it_ids, iv_ids, snp_file, f"{gwas_prefix}.res",
                                method=method_full, libpath=args.libpath
                            )
                            accuracies.append(safe_accuracy(metrics))

                        return float(np.mean(accuracies)) if accuracies else 0.0

                    best_full_params = None
                    best_full_value = None
                    for rank, t in enumerate(candidates, start=1):
                        try:
                            v = _full_eval_params(dict(t.params), int(t.number), int(rank))
                        except Exception as e:
                            LOGGER.log(f"Two-stage full-metric reevaluation failed for trial {t.number}: {e}")
                            continue
                        LOGGER.log(f"Two-stage full metric: trial {t.number} mean_acc={v:.6f} params={t.params}")
                        if best_full_value is None or v > best_full_value:
                            best_full_value = v
                            best_full_params = dict(t.params)

                    if best_full_params:
                        best_params = best_full_params
                        fold_record["two_stage"] = True
                        fold_record["two_stage_topk"] = int(len(candidates))
                        fold_record["best_params_proxy"] = dict(best_params_proxy)
                        fold_record["best_params_full"] = dict(best_params)
                except Exception as e:
                    LOGGER.log(f"Two-stage optimization failed (fallback to proxy best): {e}")

            LOGGER.log(f"Best Params (outer {outer_i}): {best_params}")
            fold_record["best_params"] = dict(best_params) if best_params else {}
            selection_method_best = selection_method_fold
            if isinstance(best_params, dict) and isinstance(best_params.get("selection_method"), str):
                if best_params["selection_method"] in {"clump", "prune", "hybrid", "clump_cached"}:
                    selection_method_best = best_params["selection_method"]
            fold_record["selection_method"] = selection_method_best
            
            # Final GWAS + Evaluate on outer test
            final_dir = os.path.join(outer_dir, "final")
            try:
                os.makedirs(final_dir, exist_ok=True)
            except Exception:
                pass
            final_prefix = os.path.join(final_dir, run_id)
            p = best_params["p_thresh"]; r2 = best_params["clump_r2"]; kb = best_params["clump_kb"]
            final_extract = None
            if args.info_file:
                try:
                    info_thr = float(best_params.get("info_thresh", 0.0))
                except Exception:
                    info_thr = 0.0
                if info_thr > 0:
                    final_extract = ensure_info_extract(args.info_file, info_thr, args.out, target_bfile=train_bfile)
            if external_qassoc:
                LOGGER.log("Final selection using external GWAS.")
                snps = select_snps_from_gwas(
                    bfile=train_bfile,
                    gwas_file=external_qassoc,
                    keep_ids_file=ot_ids,
                    out_prefix=final_prefix,
                    p_thresh=p,
                    ld_r2=r2,
                    window_kb=kb,
                    plink_threads=PLINK_THREADS,
                    extract_file=final_extract,
                    is_fastgwa=False,
                    selection_method=selection_method_best,
                    max_snps=max_snps_fold,
                )
            else:
                gwas_out_prefix = f"{final_prefix}.gwas"
                gwas_result_file, is_plink2, _ = run_gwas_only(
                    train_bfile, train_pheno, ot_ids, gwas_out_prefix, PLINK_THREADS, use_fastgwa=True, background=False
                )
                snps = select_snps_from_gwas(
                    bfile=train_bfile,
                    gwas_file=gwas_result_file,
                    keep_ids_file=ot_ids,
                    out_prefix=final_prefix,
                    p_thresh=p,
                    ld_r2=r2,
                    window_kb=kb,
                    plink_threads=PLINK_THREADS,
                    extract_file=final_extract,
                    is_fastgwa=is_plink2,
                    selection_method=selection_method_best,
                    max_snps=max_snps_fold,
                )
            LOGGER.log(f"Selected {len(snps)} SNPs.")
            fold_record["selected_snps"] = int(len(snps))
            if len(snps) > 0:
                snp_file = f"{final_prefix}.snps"
                pd.DataFrame(snps, columns=["SNP"]).to_csv(snp_file, index=False, header=False)
                final_method = resolve_method(args.method, best_params, args.proxy_metric)
                final_metrics = evaluate_performance(
                    train_bfile, train_bfile, train_pheno, train_pheno,
                    ot_ids, ov_ids, snp_file, final_prefix,
                    method=final_method, libpath=args.libpath
                )
                final_acc_raw = final_metrics.get("accuracy", float("nan"))
                final_acc = safe_accuracy(final_metrics)
                outer_accuracies.append(final_acc)
                fold_record["final_acc"] = final_acc_raw if not (isinstance(final_acc_raw, float) and np.isnan(final_acc_raw)) else None
                fold_record["final_r2"] = final_metrics.get("r2", None)
                fold_record["final_spearman"] = final_metrics.get("spearman", None)
                if isinstance(final_acc_raw, float) and np.isnan(final_acc_raw):
                    LOGGER.log(f"Outer Fold {outer_i} Test Accuracy: NA")
                else:
                    LOGGER.log(f"Outer Fold {outer_i} Test Accuracy: {final_acc}")
            else:
                outer_accuracies.append(0.0)
                fold_record["final_acc"] = 0.0
                LOGGER.log("No SNPs selected. Accuracy: 0.0")
            if fold_by_no is not None:
                fold_by_no[int(outer_i)] = fold_record
            else:
                fold_records.append(fold_record)
        
        if outer_accuracies or (fold_by_no is not None and fold_by_no):
            if fold_by_no is not None:
                fold_records = []
                for i in range(1, int(args.folds) + 1):
                    fr = fold_by_no.get(int(i))
                    if not isinstance(fr, dict):
                        fr = {"outer_fold": int(i)}
                    fold_records.append(fr)
            try:
                for t in list(BASELINE_THREADS):
                    try:
                        t.join()
                    except Exception:
                        pass
            except Exception:
                pass

            try:
                with BASELINE_LOCK:
                    baseline_snapshot = dict(BASELINE_RESULTS)
            except Exception:
                baseline_snapshot = {}

            baseline_accuracies = []
            for fr in fold_records:
                try:
                    fold_no = int(fr.get("outer_fold"))
                except Exception:
                    continue
                m_full = baseline_snapshot.get(f"outer{fold_no}:baseline_full_input")
                if isinstance(m_full, dict):
                    fr["baseline_full_acc"] = safe_accuracy(m_full)
                    fr["baseline_full_r2"] = m_full.get("r2", None)
                    fr["baseline_full_spearman"] = m_full.get("spearman", None)
                if fr.get("baseline_full_acc") is not None:
                    try:
                        baseline_accuracies.append(float(fr["baseline_full_acc"]))
                    except Exception:
                        pass

                m_orig = baseline_snapshot.get(f"outer{fold_no}:baseline_orig")
                if isinstance(m_orig, dict):
                    fr["baseline_orig_acc"] = safe_accuracy(m_orig)

                m_custom = baseline_snapshot.get(f"outer{fold_no}:baseline_custom")
                if isinstance(m_custom, dict):
                    fr["baseline_custom_acc"] = safe_accuracy(m_custom)

                m_cps = baseline_snapshot.get(f"outer{fold_no}:baseline_custom_plus_selected")
                if isinstance(m_cps, dict):
                    fr["baseline_custom_plus_selected_acc"] = safe_accuracy(m_cps)

                comm = {}
                for k, v in baseline_snapshot.items():
                    if not k.startswith(f"outer{fold_no}:baseline_comm_"):
                        continue
                    if isinstance(v, dict):
                        comm[k.split(":", 1)[1]] = safe_accuracy(v)
                if comm:
                    fr["baseline_commercial"] = comm
                    if fr.get("baseline_commercial_acc") is None:
                        try:
                            fr["baseline_commercial_acc"] = float(max(comm.values()))
                        except Exception:
                            pass

                adds = {}
                for k, v in baseline_snapshot.items():
                    if not k.startswith(f"outer{fold_no}:baseline_add_"):
                        continue
                    if isinstance(v, dict):
                        adds[k.split(":", 1)[1]] = safe_accuracy(v)
                if adds:
                    fr["baseline_additional"] = adds

            accuracies_to_write = []
            for fr in fold_records:
                v = fr.get("final_acc")
                if isinstance(v, (int, float)) and not (isinstance(v, float) and np.isnan(v)):
                    accuracies_to_write.append(float(v))
                else:
                    accuracies_to_write.append(None)

            completed_acc = [a for a in accuracies_to_write if isinstance(a, (int, float))]
            mean_outer = float(np.mean(completed_acc)) if completed_acc else 0.0
            if len(completed_acc) == int(args.folds):
                LOGGER.log(f"Mean Test Accuracy across {args.folds} folds: {mean_outer}")
            else:
                LOGGER.log(f"Mean Test Accuracy across {len(completed_acc)}/{args.folds} completed folds: {mean_outer}")

            try:
                os.makedirs(summary_dir, exist_ok=True)
            except Exception:
                pass
            summary_txt = os.path.join(summary_dir, "cv_summary.txt")
            trait = infer_trait_name_from_pheno(args.pheno or getattr(args, "train_pheno", None))
            with open(summary_txt, "w") as f:
                f.write(f"Trait: {trait}\n")
                f.write(f"RunID: {run_id}\n")
                f.write(f"OutDir: {run_dir}\n")
                f.write(f"Accuracies: {json.dumps(accuracies_to_write)}\n")
                f.write(f"Mean: {mean_outer}\n")
                if baseline_accuracies:
                    mean_base = float(np.mean(baseline_accuracies))
                    LOGGER.log(f"Mean Baseline Accuracy across {len(baseline_accuracies)} folds: {mean_base}")
                    f.write(f"Baseline_Accuracies: {json.dumps(baseline_accuracies)}\n")
                    f.write(f"Baseline_Mean: {mean_base}\n")
                f.write(f"Folds: {json.dumps(fold_records, ensure_ascii=False)}\n")
            try:
                with open(summary_json, "w") as f:
                    json.dump(
                        {
                            "trait": trait,
                            "run_id": run_id,
                            "out_dir": run_dir,
                            "outer_folds": int(args.folds),
                            "inner_folds": int(args.inner_cv_folds),
                            "accuracies": accuracies_to_write,
                            "mean_accuracy": mean_outer,
                            "baseline_accuracies": baseline_accuracies,
                            "folds": fold_records,
                        },
                        f,
                        ensure_ascii=False,
                        indent=2,
                    )
            except Exception:
                pass
        return
    
    # 1.5 Random 50k Baseline
    if not args.skip_hiblup_baseline:
        LOGGER.log("---------------------------------------------------------")
        LOGGER.log("Running Random 50k Baseline...")
        baseline_dir = os.path.join(os.path.dirname(args.out), "baselines")
        try:
            os.makedirs(baseline_dir, exist_ok=True)
        except Exception:
            pass
        try:
            launch_baseline(
                run_random_50k_baseline,
                (
                    train_bfile, test_bfile,
                    master_train_ids, master_test_ids,
                    train_pheno, test_pheno,
                    os.path.join(baseline_dir, "baseline_random50k")
                ),
                {"info_file": args.info_file, "info_thresh": args.info_thresh},
                "Random 50k Baseline"
            )
        except Exception as e:
            LOGGER.log(f"Random 50k Baseline failed to launch: {e}")
        LOGGER.log("---------------------------------------------------------")

    # ---------------------------------------------------------
    # Baselines (Standard Mode)
    # ---------------------------------------------------------

    # 1. Existing "Original Density" Baseline
    baseline_train = args.baseline_train_bfile if args.baseline_train_bfile else args.baseline_bfile
    baseline_test = args.baseline_test_bfile if args.baseline_test_bfile else (args.baseline_bfile if args.baseline_bfile else args.baseline_train_bfile)

    if not args.skip_hiblup_baseline and baseline_train and baseline_test and os.path.exists(f"{baseline_train}.bed") and os.path.exists(f"{baseline_test}.bed"):
        try:
            baseline_dir = os.path.join(os.path.dirname(args.out), "baselines")
            try:
                os.makedirs(baseline_dir, exist_ok=True)
            except Exception:
                pass
            launch_baseline(
                run_fixed_snp_baseline,
                (
                    baseline_train, baseline_test, 
                    master_train_ids, master_test_ids, 
                    train_pheno, test_pheno,
                    os.path.join(baseline_dir, "baseline_orig"),
                    None
                ),
                None,
                "Baseline Original Density (GBLUP)"
            )
        except Exception as e:
            LOGGER.log(f"Baseline GBLUP failed: {e}")

    # 2. Full-Genome GBLUP (Input Data)
    if args.baseline_full:
        try:
            baseline_dir = os.path.join(os.path.dirname(args.out), "baselines")
            try:
                os.makedirs(baseline_dir, exist_ok=True)
            except Exception:
                pass
            launch_baseline(
                run_fixed_snp_baseline,
                (
                    train_bfile, test_bfile, 
                    master_train_ids, master_test_ids, 
                    train_pheno, test_pheno,
                    os.path.join(baseline_dir, "baseline_full_input"),
                    None
                ),
                None,
                "Baseline Full-Genome GBLUP (Input Data)"
            )
        except Exception as e:
            LOGGER.log(f"Full-Genome Baseline failed: {e}")

    # 3. Custom Low-Density Chip (No Imputation)
    if args.baseline_custom_bfile or args.baseline_custom_snps:
        c_train = args.baseline_custom_bfile if args.baseline_custom_bfile else train_bfile
        c_test = args.baseline_custom_bfile if args.baseline_custom_bfile else test_bfile
        try:
            baseline_dir = os.path.join(os.path.dirname(args.out), "baselines")
            try:
                os.makedirs(baseline_dir, exist_ok=True)
            except Exception:
                pass
            launch_baseline(
                run_fixed_snp_baseline,
                (
                    c_train, c_test,
                    master_train_ids, master_test_ids,
                    train_pheno, test_pheno,
                    os.path.join(baseline_dir, "baseline_custom"),
                    args.baseline_custom_snps
                ),
                None,
                "Custom Low-Density Chip"
            )
        except Exception as e:
            LOGGER.log(f"Custom Chip Baseline failed: {e}")

    # 4. Commercial Chips
    if args.baseline_commercial_snps:
        c_train = args.baseline_commercial_bfile if args.baseline_commercial_bfile else train_bfile
        c_test = args.baseline_commercial_bfile if args.baseline_commercial_bfile else test_bfile
        
        chip_lists = args.baseline_commercial_snps.split(",")
        for idx, chip_list in enumerate(chip_lists):
            chip_list = chip_list.strip()
            if not chip_list: continue
            try:
                baseline_dir = os.path.join(os.path.dirname(args.out), "baselines")
                try:
                    os.makedirs(baseline_dir, exist_ok=True)
                except Exception:
                    pass
                run_fixed_snp_baseline(
                    c_train, c_test, 
                    master_train_ids, master_test_ids, 
                    train_pheno, test_pheno,
                    os.path.join(baseline_dir, f"baseline_comm_{idx+1}"),
                    chip_list,
                    description=f"Commercial Chip {idx+1} ({os.path.basename(chip_list)})"
                )
            except Exception as e:
                LOGGER.log(f"Commercial Chip Baseline {idx+1} failed: {e}")

    if args.fast_mode:
        p_space = [1e-5, 1e-4, 1e-3, 0.01, 0.05]
        r2_space = [0.1, 0.2, 0.5]
        kb_space = [250, 500]
    else:
        p_space = [1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
        r2_space = [0.1, 0.2, 0.5, 0.8]
        kb_space = [250, 500, 1000]
    max_snps_current = int(args.max_snps) if args.max_snps else None
    selection_method_current = args.selection_method

    need_background = bool(args.auto_search_space) or (getattr(args, "architecture_preset", "auto") != "auto") or (args.selection_method == "auto") or bool(prior)
    if need_background:
        prior_preset = None
        if prior.get("preset") in {"polygenic", "oligogenic", "hybrid"}:
            prior_preset = prior["preset"]
        elif prior.get("architecture_preset") in {"polygenic", "oligogenic", "hybrid"}:
            prior_preset = prior["architecture_preset"]

        preset = args.architecture_preset if args.architecture_preset != "auto" else prior_preset
        diag = None
        if not preset:
            diag_file = None
            diag_is_fastgwa = False
            if external_qassoc and os.path.exists(external_qassoc):
                diag_file = external_qassoc
                diag_is_fastgwa = False
            else:
                run_dir = os.path.dirname(args.out)
                run_id = os.path.basename(args.out)
                try:
                    os.makedirs(os.path.join(run_dir, "diag"), exist_ok=True)
                except Exception:
                    pass
                diag_prefix = os.path.join(run_dir, "diag", run_id)
                diag_file, diag_is_fastgwa, _ = run_gwas_only(
                    train_bfile, train_pheno, master_train_ids, diag_prefix,
                    args.plink_threads, use_fastgwa=True, background=False, reuse=args.reuse_gwas
                )
            if diag_file and os.path.exists(diag_file):
                diag = gwas_diagnostics(diag_file, is_fastgwa=diag_is_fastgwa)
                preset = infer_architecture_preset(diag)
            else:
                preset = "hybrid"

        space = biased_search_space(
            preset=preset,
            fast_mode=bool(args.fast_mode),
            base_max_snps=int(args.max_snps),
            base_clump_cap=int(args.clump_max_summary_snps),
        )
        space = apply_prior_to_space(space, prior)
        p_space = space["p_space"]
        r2_space = space["r2_space"]
        kb_space = space["kb_space"]
        max_snps_current = space["max_snps"]
        try:
            CLUMP_MAX_SUMMARY_SNPS = int(space["clump_cap"])
        except Exception:
            pass

        if args.selection_method == "auto":
            selection_method_current = choose_selection_method(space["preset"], prior=prior)

        if LOGGER:
            if diag and diag.get("n"):
                LOGGER.log(
                    f"Architecture preset: {space['preset']} (N={diag['n']}, p<1e-4={diag['counts'].get(1e-4,0)}, lambdaGC={diag.get('lambda_gc')}); selection_method={selection_method_current}"
                )
            else:
                LOGGER.log(f"Architecture preset: {space['preset']}; selection_method={selection_method_current}")
    elif args.selection_method == "auto":
        selection_method_current = "hybrid"

    selection_method_candidates_current = None
    if args.optimize_selection_method and args.selection_method == "auto":
        preset_for_mode = None
        try:
            preset_for_mode = (space.get("preset") if isinstance(space, dict) else None) if need_background else None
        except Exception:
            preset_for_mode = None
        if not preset_for_mode:
            try:
                preset_for_mode = (args.architecture_preset if getattr(args, "architecture_preset", "auto") != "auto" else None) or "hybrid"
            except Exception:
                preset_for_mode = "hybrid"
        preset_for_mode = str(preset_for_mode).lower()
        if preset_for_mode == "oligogenic":
            selection_method_candidates_current = ["clump", "hybrid"]
        elif preset_for_mode == "polygenic":
            selection_method_candidates_current = ["prune", "hybrid"]
        else:
            selection_method_candidates_current = ["clump", "prune", "hybrid"]

    # Optuna Objective
    def objective(trial):
        # Hyperparameters
        if args.method == "auto":
            if args.proxy_metric:
                method = "rrblup"
            else:
                method = trial.suggest_categorical("method", ["rrblup", "adaptive", "lightgbm"])
        else:
            method = args.method

        p_thresh = trial.suggest_categorical("p_thresh", p_space)
        clump_r2 = trial.suggest_categorical("clump_r2", r2_space)
        clump_kb = trial.suggest_categorical("clump_kb", kb_space)
        info_thresh = trial.suggest_categorical("info_thresh", [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]) if args.info_file else 0.0
        
        # Internal CV on Training Set
        n_folds_cv = args.inner_cv_folds
        kf = KFold(n_splits=n_folds_cv, shuffle=True, random_state=42)
        accuracies = []
        
        # We split 'valid_train' indices
        X_indices = np.arange(len(valid_train))
        
        for i, (tr_idx, val_idx) in enumerate(kf.split(X_indices)):
            LOGGER.log(f"[Trial {trial.number}] --- Internal CV Fold {i+1}/{n_folds_cv} ---")
            
            inner_train = valid_train.iloc[tr_idx]
            inner_val = valid_train.iloc[val_idx]
            
            # Write ID files
            it_ids = f"{args.out}.trial{trial.number}.fold{i}.train.ids"
            iv_ids = f"{args.out}.trial{trial.number}.fold{i}.val.ids"
            inner_train[["FID", "IID"]].to_csv(it_ids, sep="\t", index=False, header=False)
            inner_val[["FID", "IID"]].to_csv(iv_ids, sep="\t", index=False, header=False)
            
            # Prepare Extract File (Info Score)
            extract_file = ensure_info_extract(args.info_file, info_thresh, args.out, target_bfile=train_bfile) if (args.info_file and info_thresh > 0) else None

            # Run GWAS on Inner Train
            gwas_prefix = f"{args.out}.trial{trial.number}.fold{i}"
            if external_qassoc:
                LOGGER.log(f"Selecting using external GWAS for fold {i+1} (p={p_thresh}, r2={clump_r2}, kb={clump_kb}, info={info_thresh})")
                gwas_result_file = external_qassoc
                is_fastgwa_out = False
            else:
                LOGGER.log(f"Running GWAS on {len(inner_train)} samples...")
                gwas_out_prefix = f"{gwas_prefix}.gwas"
                gwas_result_file, is_plink2, _ = run_gwas_only(
                    train_bfile, train_pheno, it_ids, gwas_out_prefix, PLINK_THREADS, use_fastgwa=True, background=False
                )
                is_fastgwa_out = is_plink2

            selection_method_trial = selection_method_current
            if selection_method_candidates_current:
                selection_method_trial = trial.suggest_categorical("selection_method", selection_method_candidates_current)
            snps = select_snps_from_gwas(
                bfile=train_bfile,
                gwas_file=gwas_result_file,
                keep_ids_file=it_ids,
                out_prefix=gwas_prefix,
                p_thresh=p_thresh,
                ld_r2=clump_r2,
                window_kb=clump_kb,
                plink_threads=PLINK_THREADS,
                extract_file=extract_file,
                is_fastgwa=is_fastgwa_out,
                selection_method=selection_method_trial,
                max_snps=max_snps_current,
            )
            LOGGER.log(f"Selected {len(snps)} SNPs.")

            extract_file = extract_file
            
            if len(snps) == 0:
                LOGGER.log("No SNPs selected, skipping training.")
                accuracies.append(0.0)
                running_mean = float(np.mean(accuracies)) if accuracies else 0.0
                trial.report(running_mean, step=i)
                if trial.should_prune():
                    LOGGER.log(f"Pruned Trial {trial.number} at inner fold {i+1} (mean_acc={running_mean:.4f}).")
                    raise optuna.TrialPruned()
                continue
                
            snp_file = f"{gwas_prefix}.snps"
            pd.DataFrame(snps, columns=["SNP"]).to_csv(snp_file, index=False, header=False)
            
            # Train/Eval on Inner Val
            # Note: For internal CV, both train and test come from 'train_bfile'
            # because we split the training set itself.
            
            if args.proxy_metric:
                # Use Proxy (PRS or SBLUP)
                qassoc_for_proxy = external_qassoc if external_qassoc else gwas_result_file
                
                if args.proxy_method == "sblup":
                     metrics = calculate_sblup_proxy(
                        train_bfile, train_bfile,
                        it_ids, iv_ids,
                        train_pheno,
                        snp_file,
                        qassoc_for_proxy,
                        f"{gwas_prefix}.proxy",
                        gcta_path=args.gcta_path
                    )
                     LOGGER.log(f"Fold {i+1} Proxy Accuracy (SBLUP): {metrics['accuracy']:.4f}")
                else:
                    # Default: PRS
                    metrics = calculate_prs_proxy(
                        train_bfile, train_bfile,
                        it_ids, iv_ids,
                        train_pheno,
                        snp_file,
                        qassoc_for_proxy,
                        f"{gwas_prefix}.proxy"
                    )
                    LOGGER.log(f"Fold {i+1} Proxy Accuracy (PRS): {metrics['accuracy']:.4f}")
                    
                acc = metrics["accuracy"]
            else:
                # Full GBLUP/Model
                LOGGER.log("Training model...")
                metrics = evaluate_performance(
                    train_bfile, train_bfile, train_pheno, train_pheno,
                    it_ids, iv_ids, snp_file, f"{gwas_prefix}.res",
                    method=method, libpath=args.libpath
                )
                acc_raw = metrics.get("accuracy", float("nan"))
                acc = safe_accuracy(metrics)
                if isinstance(acc_raw, float) and np.isnan(acc_raw):
                    LOGGER.log(f"Fold {i+1} Accuracy: NA (Method: {method})")
                else:
                    LOGGER.log(f"Fold {i+1} Accuracy: {acc:.4f} (Method: {method})")
            
            accuracies.append(acc)
            running_mean = float(np.mean(accuracies)) if accuracies else 0.0
            trial.report(running_mean, step=i)
            if trial.should_prune():
                LOGGER.log(f"Pruned Trial {trial.number} at inner fold {i+1} (mean_acc={running_mean:.4f}).")
                raise optuna.TrialPruned()
            
            # Cleanup
            # ...
            
        mean_acc = float(np.mean(accuracies)) if len(accuracies) else 0.0
        LOGGER.log(f"Trial {trial.number} Mean Accuracy: {mean_acc:.4f}")
        return mean_acc

    # Optimization
    pruner = optuna.pruners.MedianPruner(
        n_startup_trials=max(0, int(args.optuna_startup_trials)),
        n_warmup_steps=max(1, int(args.inner_cv_folds) // 2),
    ) if args.proxy_metric else None
    callbacks = [make_optuna_early_stop_callback(args.optuna_patience, args.optuna_min_delta)]
    study = optuna.create_study(
        direction="maximize",
        sampler=optuna.samplers.TPESampler(seed=args.optuna_seed),
        pruner=pruner,
    )
    study.optimize(objective, n_trials=args.trials, n_jobs=args.n_jobs, callbacks=callbacks)
    
    best_params = study.best_params
    best_params_proxy = dict(best_params) if best_params else {}

    if args.proxy_metric and args.two_stage_opt:
        try:
            try:
                topk = int(getattr(args, "two_stage_topk", 8))
            except Exception:
                topk = 8
            topk = max(1, topk)
            completed = [
                t for t in study.trials
                if getattr(t, "value", None) is not None and getattr(t, "state", None) == optuna.trial.TrialState.COMPLETE
            ]
            completed.sort(key=lambda t: float(t.value), reverse=True)
            candidates = completed[:topk]
            if candidates:
                LOGGER.log(f"Two-stage: reevaluating top {len(candidates)} proxy trials with full metric...")

            def _full_eval_params(params: dict, trial_no: int, rank: int) -> float:
                p_thresh = params.get("p_thresh")
                clump_r2 = params.get("clump_r2")
                clump_kb = params.get("clump_kb")
                try:
                    info_thresh = float(params.get("info_thresh", 0.0))
                except Exception:
                    info_thresh = 0.0
                method_full = resolve_method(args.method, params, proxy_metric=False)
                n_folds_cv = args.inner_cv_folds
                kf = KFold(n_splits=n_folds_cv, shuffle=True, random_state=42)
                accuracies = []
                X_indices = np.arange(len(valid_train))

                run_dir = os.path.dirname(args.out)
                run_id = os.path.basename(args.out)
                stage_root = os.path.join(run_dir, "two_stage", f"trial{trial_no}", f"rank{rank}")
                try:
                    os.makedirs(stage_root, exist_ok=True)
                except Exception:
                    pass

                for i, (tr_idx, val_idx) in enumerate(kf.split(X_indices)):
                    inner_train = valid_train.iloc[tr_idx]
                    inner_val = valid_train.iloc[val_idx]

                    fold_dir = os.path.join(stage_root, f"fold{i+1}")
                    try:
                        os.makedirs(fold_dir, exist_ok=True)
                    except Exception:
                        pass

                    it_ids = os.path.join(fold_dir, "train.ids")
                    iv_ids = os.path.join(fold_dir, "val.ids")
                    inner_train[["FID", "IID"]].to_csv(it_ids, sep="\t", index=False, header=False)
                    inner_val[["FID", "IID"]].to_csv(iv_ids, sep="\t", index=False, header=False)

                    extract_file = ensure_info_extract(args.info_file, info_thresh, args.out) if (args.info_file and info_thresh > 0) else None
                    gwas_prefix = os.path.join(fold_dir, run_id)
                    if external_qassoc:
                        gwas_result_file = external_qassoc
                        is_fastgwa_out = False
                    else:
                        gwas_out_prefix = f"{gwas_prefix}.gwas"
                        gwas_result_file, is_plink2, _ = run_gwas_only(
                            train_bfile, train_pheno, it_ids, gwas_out_prefix, PLINK_THREADS, use_fastgwa=True, background=False
                        )
                        is_fastgwa_out = is_plink2

                    snps = select_snps_from_gwas(
                        bfile=train_bfile,
                        gwas_file=gwas_result_file,
                        keep_ids_file=it_ids,
                        out_prefix=gwas_prefix,
                        p_thresh=p_thresh,
                        ld_r2=clump_r2,
                        window_kb=clump_kb,
                        plink_threads=PLINK_THREADS,
                        extract_file=extract_file,
                        is_fastgwa=is_fastgwa_out,
                        selection_method=(
                            str(params.get("selection_method"))
                            if isinstance(params.get("selection_method"), str) and params.get("selection_method") in {"clump", "prune", "hybrid", "clump_cached"}
                            else selection_method_current
                        ),
                        max_snps=max_snps_current,
                    )

                    if len(snps) == 0:
                        accuracies.append(0.0)
                        continue

                    snp_file = f"{gwas_prefix}.snps"
                    pd.DataFrame(snps, columns=["SNP"]).to_csv(snp_file, index=False, header=False)
                    metrics = evaluate_performance(
                        train_bfile, train_bfile, train_pheno, train_pheno,
                        it_ids, iv_ids, snp_file, f"{gwas_prefix}.res",
                        method=method_full, libpath=args.libpath
                    )
                    accuracies.append(safe_accuracy(metrics))

                return float(np.mean(accuracies)) if accuracies else 0.0

            best_full_params = None
            best_full_value = None
            for rank, t in enumerate(candidates, start=1):
                try:
                    v = _full_eval_params(dict(t.params), int(t.number), int(rank))
                except Exception as e:
                    LOGGER.log(f"Two-stage full-metric reevaluation failed for trial {t.number}: {e}")
                    continue
                LOGGER.log(f"Two-stage full metric: trial {t.number} mean_acc={v:.6f} params={t.params}")
                if best_full_value is None or v > best_full_value:
                    best_full_value = v
                    best_full_params = dict(t.params)

            if best_full_params:
                best_params = best_full_params
                LOGGER.log(f"Two-stage: selected params by full metric; proxy_best={best_params_proxy} full_best={best_params}")
        except Exception as e:
            LOGGER.log(f"Two-stage optimization failed (fallback to proxy best): {e}")

    LOGGER.log(f"Best Params: {best_params}")
    selection_method_best = selection_method_current
    if isinstance(best_params, dict) and isinstance(best_params.get("selection_method"), str):
        if best_params["selection_method"] in {"clump", "prune", "hybrid", "clump_cached"}:
            selection_method_best = best_params["selection_method"]
    
    # Final Training & Testing
    LOGGER.log("Starting Final Evaluation...")
    
    # 1. GWAS on Full Training Set
    run_dir = os.path.dirname(args.out)
    run_id = os.path.basename(args.out)
    final_dir = os.path.join(run_dir, "final")
    try:
        os.makedirs(final_dir, exist_ok=True)
    except Exception:
        pass
    final_prefix = os.path.join(final_dir, run_id)
    p = best_params["p_thresh"]
    r2 = best_params["clump_r2"]
    kb = best_params["clump_kb"]
    final_extract = None
    if args.info_file:
        try:
            info_thr = float(best_params.get("info_thresh", 0.0))
        except Exception:
            info_thr = 0.0
        if info_thr > 0:
            final_extract = ensure_info_extract(args.info_file, info_thr, args.out, target_bfile=train_bfile)

    final_gwas_file = None
    final_is_fastgwa = False
    if external_qassoc:
        final_gwas_file = external_qassoc
        final_is_fastgwa = False
    else:
        gwas_out_prefix = f"{final_prefix}.gwas"
        final_gwas_file, final_is_fastgwa, _ = run_gwas_only(
            train_bfile, train_pheno, master_train_ids, gwas_out_prefix, PLINK_THREADS, use_fastgwa=True, background=False
        )

    snps = select_snps_from_gwas(
        bfile=train_bfile,
        gwas_file=final_gwas_file,
        keep_ids_file=master_train_ids,
        out_prefix=final_prefix,
        p_thresh=p,
        ld_r2=r2,
        window_kb=kb,
        plink_threads=PLINK_THREADS,
        extract_file=final_extract,
        is_fastgwa=final_is_fastgwa,
        selection_method=selection_method_best,
        max_snps=max_snps_current,
    )
        
    # Cap by max_snps
    cap_snps = None
    try:
        cap_snps = int(max_snps_current) if max_snps_current else None
    except Exception:
        cap_snps = None
    if cap_snps is None:
        cap_snps = int(args.max_snps) if getattr(args, "max_snps", None) else None

    if cap_snps and len(snps) > cap_snps:
        snps = snps[:cap_snps]
        LOGGER.log(f"Selected {len(snps)} SNPs (capped to max_snps={cap_snps}).")
    else:
        LOGGER.log(f"Selected {len(snps)} SNPs.")

    if len(snps) > 0:
        snp_file = f"{final_prefix}.snps"
        pd.DataFrame(snps, columns=["SNP"]).to_csv(snp_file, index=False, header=False)

        final_method = resolve_method(args.method, best_params, args.proxy_metric)
        final_metrics = evaluate_performance(
            train_bfile, test_bfile, train_pheno, test_pheno,
            master_train_ids, master_test_ids, snp_file, final_prefix,
            method=final_method, libpath=args.libpath
        )
        final_acc_raw = final_metrics.get("accuracy", float("nan"))
        final_acc = safe_accuracy(final_metrics)
        if isinstance(final_acc_raw, float) and np.isnan(final_acc_raw):
            LOGGER.log("Final Test Metrics: Acc=NA, R2=NA, Spear=NA")
        else:
            LOGGER.log(f"Final Test Metrics: Acc={final_acc:.4f}, R2={final_metrics['r2']:.4f}, Spear={final_metrics['spearman']:.4f}")

        if args.proxy_metric:
            qassoc_final = external_qassoc if external_qassoc else final_gwas_file
            if args.proxy_method == "sblup":
                pm = calculate_sblup_proxy(
                    train_bfile, test_bfile,
                    master_train_ids, master_test_ids,
                    test_pheno,
                    snp_file,
                    qassoc_final,
                    f"{final_prefix}.proxy_final",
                    gcta_path=args.gcta_path
                )
                LOGGER.log(f"Final Validation Proxy (SBLUP): {pm['accuracy']:.4f}")
            else:
                pm = calculate_prs_proxy(
                    train_bfile, test_bfile,
                    master_train_ids, master_test_ids,
                    test_pheno,
                    snp_file,
                    qassoc_final,
                    f"{final_prefix}.proxy_final"
                )
                LOGGER.log(f"Final Validation Proxy (PRS): {pm['accuracy']:.4f}")

        with open(os.path.join(final_dir, "final_result.txt"), "w") as f:
            f.write(f"Accuracy: {final_acc_raw}\n")
            f.write(f"R2: {final_metrics.get('r2', float('nan'))}\n")
            f.write(f"Spearman: {final_metrics.get('spearman', float('nan'))}\n")
            f.write(f"SNPs: {len(snps)}\n")
            f.write(f"Params: {best_params}\n")
    else:
        LOGGER.log("No SNPs selected. Accuracy: NA")

if __name__ == "__main__":
    main()
