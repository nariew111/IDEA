# IDEA-design scoring weights and calibration resources

This resource separates an **unchanged historical scoring object** from a
**new, dated five-sample calibration reconstruction**. The original study ABS
membership has not been recovered. The reconstruction must not be described as
the calibration that generated the manuscript's existing panels or results.

## Fixed scoring weights

[`weights/feature_importance1.csv`](weights/feature_importance1.csv) is a
byte-for-byte copy of the publication scoring file identified in the author's
2026-07-31 analysis lock record and called by the archived panel-design scripts:

| Feature | Historical name | Fixed weight |
| --- | --- | ---: |
| Minor allele frequency | `MAF` | 1174 |
| Local LD score | `L2` | 655 |
| Conserved-region annotation | `Conserved` | 470 |
| Conservation score | `Conserved_score` | 400 |
| CDS, intron, GC, gene density and TSS-distance features | `CDS`, `Intron`, `GC`, `GD`, `DTCT` | 0 |

The original file's SHA256 is
`584bfce3791adcfa4384be6d2529635bb9a40ce4df1acf0142a34591f242fda8`.
[`weights/retained_weights.csv`](weights/retained_weights.csv) removes only
zero-weight rows for convenient use with the four-feature calibration exports.
The nonzero values are unchanged and are not normalized or retrained.

These are preserved development-stage empirical scoring weights informed by
feature modelling. The archive does not establish a recoverable, closed-form
SHAP aggregation yielding these exact numbers. They should not be relabelled as
raw SHAP values or replaced by the separate historical `feature_importance.csv`.

## New five-sample reconstruction: `reconstructed_20261113`

At the author's request on 2026-11-13, five samples were selected per population
from the existing reference-data benchmark summaries. For each population:

1. The eligible pool is the intersection of sample IDs in every retained-feature
   summary at every listed density. All 20 populations have 30 common IDs.
2. Rank IDs by SHA256 of the UTF-8 string
   `20261113<TAB>species<TAB>population<TAB>sample_id`; select the first five.
   Selection uses IDs and the fixed seed, without consulting F1 values.
3. Keep those same five samples for every feature and density in the population.
4. Compute each bin's arithmetic mean from their archived F1 values. These inputs
   were rounded to three decimal places; no unrounded precision is recovered.
5. Order bins by mean F1. A higher mean receives a larger priority rank; ties
   are broken by ascending original bin ID, matching the existing reader.

This reconstruction reuses existing imputation-evaluation outputs. It does not
perform new genotype imputation, fit new feature weights or repeat prediction
benchmarks. The other 25 archived samples are not presented as a newly validated
independent test set. Their earlier use in development and evaluation remains
part of the historical study record.

| File | Contents |
| --- | --- |
| [`manifest.json`](manifest.json) | Protocol, fixed seed, coverage, provenance and explicit validation limits |
| [`reconstructed_20261113/selected_samples.tsv`](reconstructed_20261113/selected_samples.tsv) | The 100 selected public-reference sample IDs, five per population |
| [`reconstructed_20261113/profiles.tsv`](reconstructed_20261113/profiles.tsv) | All 140 population/density profiles and their genome assemblies |
| [`reconstructed_20261113/calibration.tsv`](reconstructed_20261113/calibration.tsv) | 2,800 feature-bin means, sample counts and priorities |
| [`source_files.tsv`](source_files.tsv) | Relative paths and SHA256 hashes for the 560 historical source summaries |
| [`historical/reader_priorities.tsv`](historical/reader_priorities.tsv) | Values and ranks actually selected by the historical summary reader |
| [`checksums.tsv`](checksums.tsv) | SHA256 hashes for the generated resource files |

No individual genotype or phenotype records, individual F1 rows, UK Biobank
participant data or company pig measurements are included. The human
SNP-imputation calibration inputs are from 1000 Genomes/IGSR, distinct from the
study's UK Biobank prediction analysis. Sample IDs identify reference-resource
records and make the new selection traceable.

## Coverage and density labels

| Species | Populations | Representative chromosome | Archived density labels | Genome build |
| --- | --- | --- | --- | --- |
| Pig | D, DJ, DT, LWU, MS, QP, TC, TP, WS, Y | 6 | 100, 200, 400, 700, 1000, 1500 | Sscrofa11.1 |
| Human | AFR, AMR, EAS, EUR, SAS | 22 | 220, 400, 600, 800, 1200, 1600 | GRCh37/hg19 |
| Chicken | Glh, Leg, RJF, Tibe, YpD | 17 | 10, 20, 30, 40, 47, 60, 80, 100, 200, 300 | GRCg6a/galGal6 |

The density is an archived representative-chromosome setting, not a requested
genome-wide marker total. Selecting a different panel size does not establish
that its performance has been benchmarked. See `profiles.tsv` for valid pairs;
the exporter rejects unknown profiles instead of silently choosing another one.

## Export a profile for IDEA-design

From the cloned repository root:

```bash
python scripts/export_calibration.py \
  --species pig --population D --density 1000 \
  --output-dir results/calibration/pig_D_1000

idea-design run \
  --feature-file path/to/annotated_candidates_with_matching_bins.tsv \
  --feature-weight results/calibration/pig_D_1000/feature_importance1.csv \
  --summary-root results/calibration/pig_D_1000 \
  --target-snp-total 10000 \
  --select-only \
  --select-out results/selected_snps.txt
```

The exporter writes four `<feature>/summary.txt` files, the four unchanged
nonzero weights and a `calibration_metadata.json` identifying the reconstruction.
For other populations, change the species, population and density together.
The generated files are compatible with both `select` and `run --select-only`;
the population-folder `select` command additionally needs the directory layout
described in the [IDEA-design guide](../../docs/idea_design.md).

**Feature bins and genomic windows are different.** The summaries preserve six
numeric-feature bin IDs and two conserved-region bin IDs and calibrate their
priority ordering. They do not contain frozen numerical feature thresholds or
per-SNP bin assignments. Supply candidate annotations with matching
`MAF_BIN`, `L2_BIN`, `Conserved_BIN` and `Conserved_score_BIN` columns when
preserving a pre-existing bin definition. The existing selectors can discretize
unbinned features themselves, but their equal-width fallback must not be assumed
to reproduce historical empirical bin assignments. These resources do not
replace candidate feature preparation, genotype inputs or reference panels.

## Historical reader archive

The historical reader uses an explicit `Average` row if present; otherwise it
uses the last data row. For the 560 density-specific source files in this bundle,
there is no explicit average row. The archived reader priorities therefore
describe the last-row behavior, not a mean across 30 samples and not a verified
five-sample ABS calibration. This limitation is recorded instead of silently
replacing the inputs and attributing a new mean to the original study.

To export the recorded historical behavior explicitly, add `--historical` to
the export command. The output row is labelled `historical_effective`, and its
metadata identifies that provenance. It is not labelled as a population average.
This archive supports auditing the reader's behavior; it does not by itself
establish which historical output panel used each file.

## Rebuild from the author's source archive

```bash
python scripts/build_calibration_resources.py \
  --source-root /path/to/author_project \
  --output-dir /tmp/rebuilt_idea_design_resources
```

The source archive must contain the paths listed in `source_files.tsv` and the
original fixed-weight file. It is not supplied by downloading the software.
The builder verifies the locked weight hash, checks the retained-feature F1
values, reproduces the deterministic sample selection and regenerates aggregate
files without exporting individual evaluation rows. These scripts use Python's
standard library. The checked-in profiles can be exported without the source
archive or participant-level data.
