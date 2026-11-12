# Installation and external tools

The source package requires Python 3.9 or later. The release preparation checks
use Linux and Python 3.10; the included CI workflow runs Python 3.11 on Linux.
Python dependencies are declared in `pyproject.toml`. Legacy data readers use
pandas APIs available before version 3, so the package specifies `pandas<3`.

## Python installation

```bash
git clone https://github.com/nariew111/IDEA.git
cd IDEA
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test]"
idea-design --help
idea-gp status
idea-gp run --help
pytest -q
```

Alternatively, run `conda env create -f environment.yml` from the repository root
and then `conda activate imputchip-idea`. The Conda environment also requests
common external executables. Conda provisioning has not been validated by the
local Python-package checks.

## External executables

| Workflow | Tools and resources |
| --- | --- |
| Small feature-table example | Python dependencies; no genomic executables |
| Feature generation | PLINK, LDSC and bedtools as required by the annotations supplied |
| Phasing, imputation and accuracy | bcftools, PLINK, Java, and user-provided Beagle and conform-gt JAR files |
| Trait prediction | PLINK/PLINK2 and the external model backend selected for the run; R and model packages, HIBLUP or GCTA may be required by the selected path |

For IDEA-GP, executable paths can be set using `IDEA_GP_PLINK`, `IDEA_GP_PLINK2`,
`IDEA_GP_RSCRIPT`, `IDEA_GP_HIBLUP` and `IDEA_GP_GCTA`, or the corresponding
runner arguments shown by `idea-gp run --help`.

The `configs/idea.example.yaml` file is accepted by `idea-design run --config`.
Its genomic inputs and JAR paths must be set for the actual analysis. IDEA-GP
uses command-line arguments; `configs/idea_gp.example.sh` is a runnable argument
template. The GP runner does not accept an IDEA-GP YAML configuration file.

## Verification scope

The included small example checks SNP selection using six artificial candidate
rows and produces three SNP IDs. Installation and that example do not establish
scientific reproducibility of imputation or phenotype-prediction benchmarks.
See [the reproducibility record](manuscript_reproducibility.md) for the remaining
study-specific materials.
