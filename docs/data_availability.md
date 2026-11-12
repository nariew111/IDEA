# Data resources and access

This software repository contains source code, documentation, tests and a small
artificial feature-table example. Its MIT software licence applies to the
software; access to third-party research data is governed by the data providers.

## Reference resources used in the associated study

| Resource | Use in the associated study | Provider entry point |
| --- | --- | --- |
| PHARP 4.0 | Porcine genotype reference resources | https://alphaindex.zju.edu.cn/PHARP/index.php/ |
| ChickenGTEx | Chicken genotype resources | https://ngdc.cncb.ac.cn/chickengtex/ |
| 1000 Genomes / IGSR | Human SNP-imputation benchmark resources | https://www.internationalgenome.org/data/ |
| HGSVC3, including Freeze4 resources | Human structural-variant benchmark targets | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/ |
| ImputeSV | Structural-variant imputation reference resources | https://yanglab.westlake.edu.cn/ImputeSV/resource |
| UK Biobank | Authorized human genotypes and phenotypes for chromosome-22 prediction analyses | https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/ |

Exact releases, genome assemblies, reference filtering and benchmark sample
identifiers must be taken from the study's Methods and accompanying data
records. A resource entry point does not identify every file used in an analysis.

## Controlled-access and company-derived data

UK Biobank individual-level genotype and phenotype data require authorization
through UK Biobank. Access is not granted by downloading this repository.

Company-derived porcine genotypes and measured phenotypes used in the study
are not publicly released. Inquiries about access can be directed to the
corresponding authors, Zhe Zhang (`zhe_zhang@zju.edu.cn`) or Yuchun Pan
(`panyuchun1963@aliyun.com`), for discussion with the data holder. This contact
route does not imply unrestricted redistribution or guaranteed access. The
study's Data availability statement must specify the actual access conditions.

These individual-level datasets are not distributed in this software repository.

## Example and study outputs

The artificial six-SNP example is documented in
[`examples/toy/README.md`](../examples/toy/README.md). Its values and weights
were manually specified without participant-level data.

The prepared manuscript supplementary tables and numerical source-data files
are separate publication materials. Their final archive and the fixed SNP lists
should be cited in the manuscript once the matching records have been deposited.
The current status of these files and the actual calibration parameters is
listed in the [reproducibility record](manuscript_reproducibility.md).
