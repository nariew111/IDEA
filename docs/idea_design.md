# IDEA-design: inputs, selection and imputation

[Back to the project overview](../README.md)

IDEA-design uses genomic features, population-specific feature-bin priorities
and fixed feature weights to select SNPs for an imputation-oriented assay panel.
The selector consumes prepared calibration summaries and weights. The
study-specific XGBoost/SHAP calibration work shown in the framework figure is
not fitted automatically by the small selection example.

## Input files

| Input | Format and interpretation |
| --- | --- |
| Candidate feature table | Tab-delimited, with `SNP` and numeric feature columns matching the weight file |
| Variant coordinates | Historical identifiers encode position as `chromosome_position`; the `run --select-only` loader also accepts explicit `CHR`/`POS` columns |
| Fixed feature weights | CSV with `feature,importance` (or `feature,weight`) |
| Feature-bin calibration | A directory containing one summary per feature, such as `MAF/summary.txt` |
| Marker budget | Total count for `run --target-snp-total`, or per-chromosome count for `run --target-snp-count` |

Use a consistent genome assembly, chromosome convention and SNP identifier
scheme across candidate features, reference genotypes and downstream extraction.
For the population-folder `select` command, use identifiers in the historical
`chromosome_position` format. Arbitrary rsIDs without coordinate columns cannot
provide the positions required for selection.

The [bundled feature table](../examples/toy/pig/toy/toy_feature_8.txt) illustrates
the column layout:

```text
SNP     MAF   L2    Conserved  Conserved_score
1_100   0.05  0.10  0          0.00
1_200   0.20  0.40  1          0.30
```

The `.txt` suffix is historical; these are tab-delimited tables.
The feature generator prepares MAF, LD-score and conservation-related columns;
other study annotations must be prepared and aligned to the candidate SNPs
before selection.

A feature-summary file contains columns such as `bin1`, `bin2`, `bin3` and an
`average` row with the empirical accuracy for each feature bin. Higher mean
accuracy gives a higher priority. The selector ranks these summaries; a missing
or unrecognized summary can trigger a fallback ranking. Provide the actual
calibration summaries when applying the calibrated method.

The four feature weights and summaries under `examples/toy/` are manually chosen
demonstration inputs. Actual archived fixed weights are supplied separately in
[`resources/idea_design/`](../resources/idea_design/README.md), together with
historical reader priorities and a new, dated five-sample calibration
reconstruction. The reconstruction is not the original study ABS. See the
resource guide and [reproducibility status](manuscript_reproducibility.md).

To export one checked-in reconstruction profile:

```bash
python scripts/export_calibration.py \
  --species human --population EUR --density 1600 \
  --output-dir results/calibration/human_EUR_1600
```

Then supply that directory's `feature_importance1.csv` and the directory itself
as `--feature-weight` and `--summary-root`. The density identifies the archived
representative-chromosome setting; it is independent of a new requested marker
total. Unknown population/density pairs are rejected.

These resources calibrate priorities for existing feature-bin labels. They do
not provide numerical feature thresholds. To retain a given bin definition,
include matching `<feature>_BIN` columns in the candidate table. Automatic
equal-width discretization of unbinned inputs must not be assumed to recreate
the study's historical empirical bin assignments.

## Select a requested marker count

For prepared tables, the `run` entry point can stop after selection:

```bash
idea-design run \
  --feature-file path/to/population_features.tsv \
  --feature-weight path/to/feature_importance.csv \
  --summary-root path/to/population_calibration \
  --target-snp-total 10000 \
  --select-only \
  --select-out results/design/selected_snps.txt
```

The output is a headerless SNP-ID list with one selected marker per line.
`--target-snp-total` requests a count across the selected chromosomes and
overrides the per-chromosome `--target-snp-count` setting. The attainable number
depends on the candidate variants and populated selection windows. Check the
output count for the actual inputs.

To restrict the candidate chromosomes, add a supported `--chr` value, for
example `--chr 1` or `--chr 1-18`. Without it, the feature table determines the
chromosomes considered. `--snp-bin-size` controls physical selection-window size
and is distinct from feature-discretization bins.

## Population-folder selection example

The smaller `select` command uses this directory layout:

```text
data_root/
  species/
    population/
      population_feature_8.txt
      feature_importance1.csv
      MAF/summary.txt
      L2/summary.txt
      Conserved/summary.txt
      Conserved_score/summary.txt
```

For the bundled artificial population:

```bash
idea-design select \
  --species pig \
  --breeds toy \
  --num-bins 3 \
  --data-root examples/toy \
  --output-dir results/toy
```

Expected output: `results/toy/toy_selected_snps.txt` containing `1_200`, `1_300`
and `1_500`, one per line. In this entry point, `--num-bins` sets the number of
position bins; it is not a universal guarantee of a genome-wide panel count.
To process multiple population folders, pass their names as one quoted string,
for example `--breeds "population_a population_b"`.

## Generate features and run imputation

Inspect feature-generation inputs before running on your own data:

```bash
idea-design generate-features --help
idea-design accuracy --help
```

The full selection-and-imputation workflow is:

```bash
idea-design run --config configs/idea.example.yaml
```

First replace the YAML's genotype, reference, benchmark, output and JAR paths.
The reference and benchmark inputs must use compatible builds and alleles, and
benchmark individuals must be separated from the reference as required by the
evaluation design. Executable and Java requirements are listed in the
[installation guide](installation.md).

The template combines artificial example features with placeholder real-data
paths to demonstrate configuration syntax. It does not include the reference
genotypes needed for an imputation run. `--skip-imputation` has a different
meaning from `--select-only`; use the latter when you only want a SNP list.

The [web resource](https://alphaindex.zju.edu.cn/IDEA) is an alternative interface
for the IDEA-design panel-generation component.
