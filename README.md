<img width="992" height="1057" alt="Gemini_Generated_Image_umz1fhumz1fhumz1" src="https://github.com/user-attachments/assets/745bd8b8-134b-4752-ad87-c1107bdc3fd1" /># IDEA: A Future-Proof Imputation-Driven SNP Array Design and Genomic Prediction Framework

IDEA is a future-proof, imputation-driven framework for SNP array design and downstream genomic prediction. It integrates two independent suites under a single, reproducible command-line entrypoint.

## Components

- ULDSuite: low-density SNP selection and imputation benchmarking (`idea/ULDSuite/2_ULDSuite_acc.py`)
- HAGSuite: genomic prediction optimization and evaluation (`idea/HAGSuite/optimize_gs.py`)

## Requirements

- Python >= 3.9
- Common external tools: plink2, bcftools
- ULDSuite external dependencies: Java (openjdk) and user-provided Beagle / conform-gt JARs (passed via CLI)
- HAGSuite optional dependencies: lightgbm, pandas-plink, matplotlib

## Installation

### Option A: Conda (recommended)

```bash
conda env create -f environment.yml
conda activate idea
pip install -e .
```
<img width="992" height="1057" alt="Gemini_Generated_Image_umz1fhumz1fhumz1" src="https://github.com/user-attachments/assets/fed9300e-5aee-4625-adff-19d9c5860f3e" />

### Option B: pip (from source)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
pip install -e ".[hag]"
```

Notes:

- `pip install -e ".[hag]"` installs common optional dependencies for HAGSuite
- External tools (plink2/bcftools/openjdk) are not Python dependencies; install them via conda or your system package manager

## Quick Start

```bash
idea --help
idea uld -- -h
idea hag -- -h
```

## Command Structure

IDEA forwards all arguments after `--` to the underlying suite script:

```bash
idea <uld|hag|github> -- [arguments forwarded to the suite...]
```

## ULDSuite (SNP Selection + Imputation Benchmark)

### What it does

1. Loads feature definitions and feature importance weights
2. Builds bin priorities using per-feature summary files (under `--summary-root`)
3. Selects a target number of SNPs per chromosome
4. Optionally runs chunked Beagle imputation and evaluates accuracy against a benchmark

### Inputs

- Feature file: `--feature-file`
- Feature weight CSV: `--feature-weight`
- Feature summary root directory: `--summary-root`
- Reference VCF: `--reference-vcf`
- Benchmark: `--benchmark`
- Tools: `--plink`, `--bcftools`, `--beagle-jar`, `--conform-jar`
- Chromosomes: `--chr` (comma-separated, or a range-style single value depending on your setup)

### Recommended usage (via YAML config)

```bash
idea uld -- --config path/to/config.yaml
```

Minimal example `config.yaml` (adjust paths to your environment):

```yaml
feature_file: path/to/features.txt
feature_weight: path/to/feature_weight.csv
summary_root: path/to/feature_summaries/
reference_vcf: path/to/reference.vcf.gz
benchmark: path/to/benchmark.vcf.gz
output_dir: ./imputation_output
beagle_jar: path/to/beagle.jar
conform_jar: path/to/conform-gt.jar
accuracy_script: path/to/accuracy_script.py
plink: plink2
bcftools: bcftools
chr: "1,2,3"
threads: 8
target_snp_count: 1000
```

### Selection-only mode

If you only want a SNP list without running imputation:

```bash
idea uld -- --config path/to/config.yaml --select-only --select-out selected_snps.txt
```

### Direct script execution (equivalent)

```bash
python idea/ULDSuite/2_ULDSuite_acc.py --config path/to/config.yaml
```

## HAGSuite (Genomic Prediction Optimization)

### What it does

HAGSuite evaluates genomic prediction performance using multiple methods (e.g., rrBLUP / lasso / KAML / LightGBM) and supports cross-validation, optional GWAS-based SNP selection, and Optuna-based hyperparameter optimization.

### Minimal run

At minimum, specify an output prefix:

```bash
idea hag -- --out results/run1
```

In practice you will typically also provide genotype and phenotype inputs. Common patterns:

### Single dataset (bfile + pheno)

```bash
idea hag -- \
  --bfile path/to/genotype_prefix \
  --pheno path/to/pheno.txt \
  --out results/run1 \
  --method auto
```

### Train/test split inputs

```bash
idea hag -- \
  --train_bfile path/to/train_prefix \
  --test_bfile path/to/test_prefix \
  --train_pheno path/to/train.pheno \
  --test_pheno path/to/test.pheno \
  --out results/run_split \
  --method rrblup
```

### Cross-validation + Optuna

```bash
idea hag -- \
  --bfile path/to/genotype_prefix \
  --pheno path/to/pheno.txt \
  --out results/cv_optuna \
  --folds 5 \
  --trials 200 \
  --method auto
```

### Direct script execution (equivalent)

```bash
python idea/HAGSuite/optimize_gs.py -h
```

## Repository Layout

```text
.
├── idea/
│   ├── ULDSuite/
│   ├── HAGSuite/
│   └── idea.py
├── environment.yml
└── pyproject.toml
```

## Publish to GitHub

Print recommended commands:

```bash
idea github --repo https://github.com/<your-username>/idea.git
```

Or run manually:

```bash
git init
git add .
git commit -m "Init IDEA: ULDSuite + HAGSuite"
git branch -M main
git remote add origin https://github.com/<your-username>/idea.git
git push -u origin main
```

## Troubleshooting

- Missing LightGBM / pandas-plink: install with `pip install -e ".[hag]"` or use the conda environment.
- Full CLI options: run `idea uld -- -h` or `idea hag -- -h`.
