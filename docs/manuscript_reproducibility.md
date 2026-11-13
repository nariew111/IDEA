# Study reproducibility record

The software package exposes the preserved IDEA-design and IDEA-GP engines.
The version in `pyproject.toml` is a package version; an exact Git commit should
be cited when identifying a software snapshot. The original local source import
is recorded by hashes in `source_snapshot.json`.

| Material | Status in this repository |
| --- | --- |
| IDEA-design feature preprocessing, selection and imputation engine | Source under `src/idea/` |
| IDEA-GP selection and prediction engine | Source under `src/idea_gp/` |
| Installation, dependencies and example commands | README, `environment.yml`, `pyproject.toml`, `configs/` |
| Detailed input/output guides and Chinese onboarding | `docs/idea_design.md`, `docs/idea_gp.md`, `docs/quickstart_zh.md` |
| Framework figure supplied with the manuscript | `assets/IDEA_framework.jpeg`, displayed in the README |
| Artificial IDEA-design selection example and expected SNPs | `examples/toy/`, `scripts/run_toy_idea_design.sh` |
| Package and CLI regression checks | `tests/`, `.github/workflows/tests.yml` |
| Fixed scoring weights called by the archived panel scripts | Original file copied unchanged to `resources/idea_design/weights/feature_importance1.csv`; checksum and provenance documented |
| Historical summary-reader priorities | `resources/idea_design/historical/reader_priorities.tsv`; actual row-selection behavior preserved |
| New five-sample calibration reconstruction | 20 populations, 140 profiles, sample IDs and aggregate bin means under `resources/idea_design/reconstructed_20261113/`; newly generated, not an original ABS record |
| Original five-sample ABS membership and numerical bin thresholds | Not recovered or verified; the new reconstruction does not resolve this historical provenance gap |
| Exact benchmark runner scripts, comparator settings, statistical summaries and figure-generation scripts | The imported software bundle maintains these separately; their final manuscript-matched collection is pending |
| Exact evaluated SNP panel lists | Final manuscript-matched files and archive location still need to be identified |
| Individual-level UK Biobank and company pig data | Access through the respective provider under the manuscript's declared conditions |
| Full IDEA-GP example with redistributable genotypes | Not included in this software snapshot |

## Interpreting the framework figure

The illustration describes the scientific workflow, including empirical
imputation calibration and feature-contribution modelling. The package's
selection commands consume prepared feature tables, bin-accuracy summaries and
feature weights. The current source bundle does not include a self-contained
script that regenerates all of the manuscript's XGBoost/SHAP calibration models
from their original empirical inputs.

The small example demonstrates the interfaces with artificial inputs. A public
example being executable does not imply that the study's parameters, genotypes
or full benchmark results are bundled with it.

## Materials needed to reproduce the manuscript

| Stage | Manuscript-matched materials needed in addition to the core code |
| --- | --- |
| Empirical calibration | Population labels, input/reference versions, feature definitions, training-panel settings, calibration summaries and fixed feature weights |
| Sparse-panel evaluation | Final SNP-ID lists, genome builds, reference/target partitioning, imputation-tool versions and comparator panel definitions |
| Trait prediction | Authorized genotype and phenotype inputs, cohort filters, fold assignments, training-only evidence, search settings, model versions and seeds |
| Figures and statistics | Numerical source data, statistical comparison definitions and plotting scripts corresponding to the final figures |

The fixed publication scoring object is now supplied as an unchanged,
checksum-identified file in the [calibration resource bundle](../resources/idea_design/README.md).
The similarly named file under `examples/toy/` remains an artificial example.
An exact study panel archive and complete figure-reproduction command are not
provided at present.

## Calibration added on 2026-11-13

The original five-sample ABS membership could not be established from the
available records. At the author's request, a new deterministic five-sample
subset was selected per population from 30 archived reference-data benchmark
samples. The same five IDs are used at all densities; selection depends on a
fixed seed and IDs, not F1 values. Archived rounded sample-level F1 values were
then averaged to create new feature-bin means and priorities. No new imputation
or prediction runs were performed, and the other 25 records are not asserted to
form a newly independent validation set.

The separate historical reader archive preserves which values the code reads
from the source files. All 560 density-specific summaries lack an explicit
`Average` row; the historical reader therefore takes their last row. Those
values are labelled as historical reader behavior, not a population mean or a
verified ABS calibration. The reconstruction does not retroactively establish
that the existing study panels or published performance used these new means.

The export script makes either resource usable with the existing selector
without changing the engine. Numerical feature-bin boundaries and annotated
candidate SNPs remain separate inputs. See the resource guide for the full
protocol, scope and usage.

For a complete manuscript reproduction record, retain input versions, genome
builds, calibration data, fixed weights, seeds, outer and inner folds, model
settings, comparator definitions and expected output files alongside the
corresponding analysis scripts. Resource entry points are listed in
[Data resources and access](data_availability.md).

The small bundled example verifies package operation and feature-based SNP
selection. Full biological benchmarks require their study-specific inputs and
external executables.
