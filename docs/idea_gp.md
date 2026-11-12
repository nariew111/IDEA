# IDEA-GP: phenotype-specific marker selection and prediction

[Back to the project overview](../README.md)

IDEA-GP operates on a candidate genotype set and selects prediction markers for
one phenotype at a time. The runner supports nested cross-validation, LD-based
selection, optional imputation-quality filtering and several prediction
backends. The command is `idea-gp`; the preserved engine may identify itself as
HAGSSuite in help text and logs.

## Required inputs

| Input | Expected format |
| --- | --- |
| Genotypes | PLINK binary prefix; supply the common prefix of `.bed`, `.bim` and `.fam`, without an extension |
| Phenotype | Whitespace-delimited file with a header and the three columns `FID`, `IID`, `Pheno` |
| Output | A writable output path supplied with `--out` |
| Model tools | PLINK/PLINK2 and the selected model's R packages or executables |

For example, the phenotype layout is:

```text
FID       IID       Pheno
SYN_F01   SYN_P01   0.20
SYN_F02   SYN_P02   -0.10
SYN_F03   SYN_P03   0.75
```

These identifiers and values are invented solely to illustrate file structure;
they are not participant records and are not a runnable genotype/phenotype
dataset. Match `FID` and `IID` to the genotype `.fam` file. Prepare one measured
trait per run, with usable numeric phenotype values. The recommended three-column
header avoids ambiguity: when multiple trait columns are supplied, some runner
paths use `Pheno`, and otherwise select the last column.

Use input genotypes and phenotypes that you are authorized to process. The
repository does not contain individual-level UK Biobank or company pig data.
See [data resources and access](data_availability.md).

## Example invocation

```bash
idea-gp run \
  --bfile path/to/genotypes \
  --pheno path/to/phenotypes.tsv \
  --out results/gp \
  --folds 5 \
  --inner_cv_folds 5 \
  --trials 20 \
  --method rrblup \
  --optuna_seed 2026
```

This is equivalent to the provided shell template:

```bash
bash configs/idea_gp.example.sh path/to/genotypes path/to/phenotypes.tsv results/gp
```

The shell file is an argument template, not a complete biological demonstration.
The runner accepts command-line flags rather than an IDEA-GP YAML file.
Inspect the options with `idea-gp run --help`.

## Cross-validation and training evidence

| Option | Meaning |
| --- | --- |
| `--folds` | Number of outer cross-validation folds |
| `--inner_cv_folds` | Number of inner folds used during optimization |
| `--trials` | Maximum number of optimization trials per outer fold |
| `--optuna_seed` | Optimization seed; record this together with all split settings |
| `--selection_method` | Selection path: `auto`, `clump`, `clump_cached`, `prune` or `hybrid` |
| `--two_stage_opt`, `--two_stage_topk` | Enable proxy screening followed by full-model verification of selected candidates |
| `--proxy_metric`, `--proxy_method` | Use a PRS or SBLUP proxy during inner optimization; the SBLUP path needs GCTA |

The method's evaluation principle is to reserve outer-test individuals for
held-out assessment, with marker-selection evidence and tuning derived from
the appropriate training partitions. Keep any externally supplied association
evidence independent of held-out individuals. Specialized options such as
`--external_assoc`, `--gwas_bfile`, `--reuse_gwas` and
`--lock_params_across_folds` change the evidence or reuse protocol; they should
not be enabled without checking the intended evaluation design.

The example above supplies generic settings. It does not specify the
manuscript's final cohort filters, fold assignments, covariate processing,
search spaces or comparator definitions.

## Optional imputation-quality scores

Use `--info_file` to supply variant-level quality information and `--info_thresh`
to set its default threshold. A simple format is:

```text
SNP       INFO
1_100     0.95
1_200     0.88
```

These rows are artificial format examples. The reader also recognizes several
identifier and quality-column aliases, including `ID`, `rsid`, `DR2` and `R2`.
Match variant identifiers to the genotype marker IDs. Quality information is
optional; its presence and definition should be recorded for each analysis.

## Prediction backends and comparison modes

The `--method` choices are `auto`, `rrblup`, `lasso`, `kaml`, `adaptive` and
`lightgbm`. Install the backend required for the selected path. R-backed methods
use the packaged `calc_gs.R` helper; other paths can require HIBLUP or GCTA.
Executable configuration is documented in [installation](installation.md).

The runner also exposes full-candidate, fixed-panel and additional baseline
options, including `--baseline_full`, `--baseline_custom_snps`,
`--baseline_commercial_snps` and `--additional_baselines`. Enabling an option
does not recreate the manuscript comparison: candidate variants, cohort splits,
model settings and comparator SNP lists must also match.

## Output locations

By default, the runner organizes files beneath the requested run root using
the phenotype-derived trait name and the run identifier. Subdirectories include
`logs/`, `folds/`, `tmp/` and, for cross-validation runs, `summary/`.
Cross-validation summaries include `cv_summary.txt` and `cv_summary.json`,
with fold results and aggregate accuracy information. Selected-SNP and
optimization outputs depend on the chosen execution path. Read the run log to
identify the resolved paths. `--no_organize_outputs` disables automatic layout.

Runtime and memory depend on sample size, candidate-marker count, fold count,
optimization trials and model backend. No end-to-end runtime is claimed for
the placeholder command. Package checks and the artificial design example are
documented separately in the [reproducibility record](manuscript_reproducibility.md).
