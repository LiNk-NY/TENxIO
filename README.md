# TENxIO

<!-- badges: start -->
<!-- badges: end -->

`TENxIO` provides methods for importing common 10X Genomic files into familiar Bioconductor classes such as `SingleCellExperiment` and `RaggedExperiment`.

## Installation

You can install the development version of TENxIO from [GitHub](https://github.com/waldronlab/TENxIO) with:

``` r
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("waldronlab/TENxIO")
```

## Supported Formats

| **Extension**       | **Class**     | **Imported as**      |
|---------------------|---------------|----------------------|
| .h5                 | TENxH5        | SingleCellExperiment |
| .mtx / .mtx.gz      | TENxMTX       | dgCMatrix            |
| .tar.gz             | TENxFileList  | SingleCellExperiment |
| peak_annotation.tsv | TENxPeaks     | GRanges              |
| fragments.tsv.gz    | TENxFragments | RaggedExperiment     |
| .tsv / .tsv.gz      | TSVFile\*     | tibble               |

Note (\*). The `TSVFile` class is used internally and not exported.

## TENxH5

``` r
library(TENxIO)
h5f <- "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
con <- TENxH5(h5f)
import(con)
```
