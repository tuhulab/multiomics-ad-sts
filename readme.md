# Multi-omics AD PhD STS (skin tape stripping)

This repo documents the  data analysis pipeline in reproducible R markdown for Multi-omics AD PhD skin tape project, submitted as [**Hu, T et al.  Investigation of the Atopic Dermatitis Epidermal Transcrip-tome by Tape Stripping and Ultra-low Cost BRB-seq . Int. J. Mol. Sci. 2022, 23, x**](https://doi.org/10.3390/xxxxx)
## Data
Data is deposited on GEO under the accession GSExxxxx (pending approval).
## Code
| Rmd                             | description                                                                                     |
|---------------------------------|-------------------------------------------------------------------------------------------------|
| 01-mapping.Rmd                  | BRB-seq data mapping                                                                            |
| 02-dataQCcuration.Rmd           | Data QC, curation and explorative analysis (RNA yield, heatmap, PCA, etc)                       |
| 03-benchmark(LSvsNL).Rmd        | Benchmark different count transformation and differential expression testing methods (LS vs Nl) |
| 04-benchmark(LSvsHC_NLvsHC).Rmd | Benchmark ... (LS vs HC, NL vs HC)                                                              |
| 05-AD_signature.Rmd             | Generate and visualize AD signatures                                                            |
| 06-BiopsyvsTape.Rmd             | Comparison between full-thickness biopsy and tape stripped skin samples                         |
| 07-data_submission_GEO.Rmd      | Submitting data to GEO                                                                          |
| 08-benchmark_table_Rpubs.Rmd    | Submitting benchmark table to Rpubs                                                             |

