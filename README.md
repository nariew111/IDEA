# IDEA: Imputation-Driven SNP Array Design and Trait-Adaptive Genomic Prediction

[![Tests](https://github.com/nariew111/IDEA/actions/workflows/tests.yml/badge.svg)](https://github.com/nariew111/IDEA/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

IDEA separates two marker-selection decisions: which SNPs to assay for genotype
recovery, and which recovered variants to use for predicting a phenotype.
**IDEA-design** selects sparse genotyping panels using genomic features and
population-specific imputation calibration. **IDEA-GP** selects trait-adaptive
prediction markers using association evidence, LD and optional imputation-quality
information within a nested cross-validation workflow.

This repository provides the two software modules, command-line entry points,
installation instructions and a small reproducible IDEA-design selection example.
The IDEA-design web resource is available at
[alphaindex.zju.edu.cn/IDEA](https://alphaindex.zju.edu.cn/IDEA).

## Framework

<img width="550" alt="IDEA framework overview" src="https://github.com/user-attachments/assets/92e24ff5-340e-4340-b5db-50866ee17104" />

The original engines are retained under `src/idea/legacy/` and
`src/idea_gp/legacy/`. Their earlier ULDSuite and HAGSSuite names may appear in
source comments and runner output. The installed command names are
`idea-design` and `idea-gp`.

## Install

```bash
git clone https://github.com/nariew111/IDEA.git
cd IDEA
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test]"
idea-design --help
idea-gp status
```

A Conda environment is also provided:

```bash
conda env create -f environment.yml
conda activate imputchip-idea
```

See [installation and external tools](docs/installation.md) for genomic
executables, JAR files and model backends required by full analyses.

## Run the small IDEA-design example

From the repository root:

```bash
bash scripts/run_toy_idea_design.sh
```

The file `results/toy/toy_selected_snps.txt` should contain:

```text
1_200
1_300
1_500
```

This artificial six-candidate example checks feature-bin ranking and SNP
selection. It does not require PLINK, a reference panel or biological data.
To choose another output location:

```bash
OUT_DIR=/tmp/idea_toy_results bash scripts/run_toy_idea_design.sh
```

## Use IDEA-design

| Command | Purpose |
| --- | --- |
| `idea-design generate-features --help` | Prepare feature tables from genomic inputs |
| `idea-design select --help` | Select from prepared population-specific feature tables |
| `idea-design run --help` | Select markers and run the imputation/evaluation workflow |
| `idea-design accuracy --help` | Evaluate recovered genotypes |
| `idea-design run-legacy --help` | Access the earlier pipeline entry point |

For the selection-and-imputation workflow, edit the input, output and tool paths
in `configs/idea.example.yaml`, then run:

```bash
idea-design run --config configs/idea.example.yaml
```

The prepared-feature selection command is used by the bundled small example.
Population calibration summaries and feature weights are inputs to the selector;
the artificial example weights are demonstration values.

## Use IDEA-GP

Inspect the underlying runner options:

```bash
idea-gp run --help
```

For your own authorized genotype and phenotype inputs:

```bash
bash configs/idea_gp.example.sh path/to/genotypes path/to/phenotypes.tsv results/gp
```

Or pass runner arguments directly:

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

`--bfile` is a PLINK binary prefix. The phenotype file must use the identifier
and phenotype format accepted by the runner. The runner consumes command-line
options; the shell example supplies them directly. Example settings are not a
locked copy of the manuscript benchmark configurations.

## Tests and reproduction materials

```bash
pytest -q
```

The tests cover command-line entry points and expected small-example selection.
The status of manuscript-specific calibration data, fixed panel lists,
benchmark scripts and external data access is documented in the
[study reproducibility record](docs/manuscript_reproducibility.md).
[Data resources and access](docs/data_availability.md) and
[Code availability wording](docs/code_availability.md) describe the scope of the
materials provided here.

## Layout

```text
src/idea/                 IDEA-design and preserved engine scripts
src/idea_gp/              IDEA-GP and preserved prediction engine
configs/                  IDEA-design YAML and IDEA-GP shell example
examples/toy/             small artificial inputs and expected SNP list
scripts/                  small-example runner
tests/                    package regression checks
docs/                     installation, data access and reproduction status
CITATION.cff              software citation metadata
LICENSE                   MIT software licence
```

## Citation and licence

Software metadata are provided in [CITATION.cff](CITATION.cff). Identify the exact
Git commit or published release when citing a software snapshot. The source
package is distributed under the [MIT License](LICENSE).
