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
| Artificial IDEA-design selection example and expected SNPs | `examples/toy/`, `scripts/run_toy_idea_design.sh` |
| Package and CLI regression checks | `tests/`, `.github/workflows/tests.yml` |
| Population-specific calibration inputs and final feature weights used in the manuscript | Study-specific files still need to be matched and deposited |
| Exact benchmark runner scripts, comparator settings, statistical summaries and figure-generation scripts | The imported software bundle maintains these separately; their final manuscript-matched collection is pending |
| Exact evaluated SNP panel lists | Final manuscript-matched files and archive location still need to be identified |
| Individual-level UK Biobank and company pig data | Access through the respective provider under the manuscript's declared conditions |
| Full IDEA-GP example with redistributable genotypes | Not included in this software snapshot |

For a complete manuscript reproduction record, retain input versions, genome
builds, calibration data, fixed weights, seeds, outer and inner folds, model
settings, comparator definitions and expected output files alongside the
corresponding analysis scripts. Resource entry points are listed in
[Data resources and access](data_availability.md).

The small bundled example verifies package operation and feature-based SNP
selection. Full biological benchmarks require their study-specific inputs and
external executables.
