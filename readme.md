# Multi-omics AD PhD project: skin tape stripping

This repo documents the  data analysis pipeline in reproducible R markdown files for Multi-omics AD PhD skin tape project, submitted as [**Hu, T et al.  Profiling the Atopic Dermatitis Epidermal Transcriptome by Tape Stripping and BRB-seq . Int. J. Mol. Sci. 2022, 23 (11), 6140**](https://doi.org/10.3390/ijms23116140)
## Data
Data is deposited on GEO under the accession [GSE199046](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199046).
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

