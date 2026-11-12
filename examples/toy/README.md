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
