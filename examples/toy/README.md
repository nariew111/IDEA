# Small IDEA-design example

This directory contains a tiny public test dataset that can run quickly on a
reviewer laptop or CI runner.

Smoke command:

```bash
scripts/run_toy_idea_design.sh
```

Expected output: `results/toy/toy_selected_snps.txt`:

```text
1_200
1_300
1_500
```

`expected_selected_snps.txt` records the same output. This small artificial
feature table checks selection and packaging; it does not reproduce biological
imputation or phenotype-prediction performance. No external genomic tools are
needed for this example.

## How the inputs were constructed

The six candidate SNP identifiers and feature values were manually specified to
make a small, deterministic selection example. Feature weights and the feature-bin
accuracy rows were also chosen by hand. They were not calculated from UK Biobank,
company pig records or other participant-level genotypes or phenotypes.

These files illustrate the expected input structure and enable installation and
selection checks. They are not empirical calibration results and must not be
used as the paper's population-specific weights or performance estimates.

| File | Purpose |
| --- | --- |
| `pig/toy/toy_feature_8.txt` | Six artificial candidate-feature rows |
| `pig/toy/feature_importance1.csv` | Four illustrative feature weights |
| `pig/toy/*/summary.txt` | Artificial feature-bin priorities |
| `expected_selected_snps.txt` | Expected output of the population-folder selector |

The marker-count `run --select-only` path and the population-folder `select`
path expose different selection settings. Their commands and input requirements
are documented in the [IDEA-design guide](../../docs/idea_design.md).
