
<!-- badges: start -->
<!-- badges: end -->

# Introduction

`TENxIO` allows users to import 10X pipeline files into known
Bioconductor classes. The package is not comprehensive, there are files
that are not supported. It currently does not support Visium datasets.
It does replace some functionality in `DropletUtils`. If you would like
a file format to be supported. Please open an issue at
<https://github.com/waldronlab/TENxIO>.

# Installation

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("waldronlab/TENxIO")
```

# Load the package

``` r
library(TENxIO)
```

# Description

`TENxIO` offers an set of classes that allow users to easily work with
files typically obtained from the 10X Genomics website. These are
commonly outputs of the Cell Ranger pipeline.

# Supported Formats

| **Extension**       | **Class**     | **Imported as**      |
|---------------------|---------------|----------------------|
| .h5                 | TENxH5        | SingleCellExperiment |
| .mtx / .mtx.gz      | TENxMTX       | dgCMatrix            |
| .tar.gz             | TENxFileList  | SingleCellExperiment |
| peak_annotation.tsv | TENxPeaks     | GRanges              |
| fragments.tsv.gz    | TENxFragments | RaggedExperiment     |
| .tsv / .tsv.gz      | TSVFile\*     | tibble               |

**Note** (\*). The `TSVFile` class is used internally and not exported.

# TENxFile

The `TENxFile` class is the catch-all class superclass that allows
transition to subclasses pertinent to specific files. It inherits from
the `BiocFile` class and allows for easy dispatching `import` methods.

``` r
showClass("TENxFile")
#> Class "TENxFile" [package "TENxIO"]
#> 
#> Slots:
#>                                                                                                                               
#> Name:                extension                  colidx                  rowidx                  remote                resource
#> Class:               character                 integer                 integer                 logical character_OR_connection
#> 
#> Extends: "BiocFile"
#> 
#> Known Subclasses: "TSVFile", "TENxFragments", "TENxH5", "TENxMTX", "TENxPeaks"
```

## `ExperimentHub` resources

`TENxFile` can handle resources from `ExperimentHub` with careful
inputs. For example, one can import a `TENxBrainData` dataset via the
appropriate `ExperimentHub` identifier (`EH1039`):

``` r
hub <- ExperimentHub::ExperimentHub()
#> snapshotDate(): 2022-10-03
hub["EH1039"]
#> ExperimentHub with 1 record
#> # snapshotDate(): 2022-10-03
#> # names(): EH1039
#> # package(): TENxBrainData
#> # $dataprovider: 10X Genomics
#> # $species: Mus musculus
#> # $rdataclass: character
#> # $rdatadateadded: 2017-10-26
#> # $title: Brain scRNA-seq data, 'HDF5-based 10X Genomics' format
#> # $description: Single-cell RNA-seq data for 1.3 million brain cells from E18 mice. 'HDF5-based 10X Genomics' format originally provided by TENx Gen...
#> # $taxonomyid: 10090
#> # $genome: mm10
#> # $sourcetype: HDF5
#> # $sourceurl: http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5
#> # $sourcesize: NA
#> # $tags: c("SequencingData", "RNASeqData", "ExpressionData", "SingleCell") 
#> # retrieve record with 'object[["EH1039"]]'
```

Currently, `ExperimentHub` resources do not have an extension and it is
best to provide that to the `TENxFile` constructor function.

``` r
fname <- hub[["EH1039"]]
TENxFile(fname, extension = "h5", group = "mm10", version = "2")
```

Note. `EH1039` is a large \~ 4GB file.

## TENxH5

`TENxIO` mainly supports version 3 and 2 type of H5 files. These are
files with specific groups and names as seen in `h5.version.map`, an
internal `data.frame` map that guides the import operations.

``` r
TENxIO:::h5.version.map
#>   Version           ID         Symbol                   Type             Ranges
#> 1       3 /features/id /features/name /features/feature_type /features/interval
#> 2       2       /genes    /gene_names                   <NA>               <NA>
```

In the case that, there is a file without genomic coordinate
information, the constructor function can take an `NA_character_` input
for the `ranges` argument.

The `TENxH5` constructor function can be used on either version of these
H5 files. In this example, we use a subset of the PBMC granulocyte H5
file obtained from the [10X
website](https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5).

``` r
h5f <- system.file(
    "extdata", "pbmc_granulocyte_ff_bc_ex.h5",
    package = "TENxIO", mustWork = TRUE
)
library(rhdf5)
h5ls(h5f)
#>               group          name       otype  dclass dim
#> 0                 /        matrix   H5I_GROUP            
#> 1           /matrix      barcodes H5I_DATASET  STRING  10
#> 2           /matrix          data H5I_DATASET INTEGER   2
#> 3           /matrix      features   H5I_GROUP            
#> 4  /matrix/features _all_tag_keys H5I_DATASET  STRING   2
#> 5  /matrix/features  feature_type H5I_DATASET  STRING  10
#> 6  /matrix/features        genome H5I_DATASET  STRING  10
#> 7  /matrix/features            id H5I_DATASET  STRING  10
#> 8  /matrix/features      interval H5I_DATASET  STRING  10
#> 9  /matrix/features          name H5I_DATASET  STRING  10
#> 10          /matrix       indices H5I_DATASET INTEGER   2
#> 11          /matrix        indptr H5I_DATASET INTEGER  11
#> 12          /matrix         shape H5I_DATASET INTEGER   2
```

Note. The `h5ls` function gives you an overview of the structure of the
file. It matches version 3 in our version map.

The show method gives an overview of the data components in the file:

``` r
con <- TENxH5(h5f)
con
#> TENxH5 object 
#> resource: /usr/local/lib/R/host-site-library/TENxIO/extdata/pbmc_granulocyte_ff_bc_ex.h5 
#> projection: SingleCellExperiment 
#> dim: 10 10 
#> rownames: ENSG00000243485 ENSG00000237613 ... ENSG00000286448 ENSG00000236601 
#> rowData names(3): ID Symbol Type 
#>   Type: Gene Expression 
#> colnames: AAACAGCCAAATATCC-1 AAACAGCCAGGAACTG-1 ... AAACCGCGTGAGGTAG-1 AAACGCGCATACCCGG-1
```

### import TENxH5 method

We can simply use the import method to convert the file representation
to a Bioconductor class representation, typically a
`SingleCellExperiment`.

``` r
import(con)
#> class: SingleCellExperiment 
#> dim: 10 10 
#> metadata(0):
#> assays(1): counts
#> rownames(10): ENSG00000243485 ENSG00000237613 ... ENSG00000286448 ENSG00000236601
#> rowData names(3): ID Symbol Type
#> colnames(10): AAACAGCCAAATATCC-1 AAACAGCCAGGAACTG-1 ... AAACCGCGTGAGGTAG-1 AAACGCGCATACCCGG-1
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: Gene Expression
#> altExpNames(0):
```

**Note**. Although the main representation in the package is
`SingleCellExperiment`, there could be a need for alternative data class
representations of the data. The `projection` field in the `TENxH5` show
method is an initial attempt to allow alternative representations.

# TENxMTX

Matrix Market formats are also supported (`.mtx` extension). These are
typically imported as SummarizedExperiment as they usually contain count
data.

``` r
mtxf <- system.file(
    "extdata", "pbmc_3k_ff_bc_ex.mtx",
    package = "TENxIO", mustWork = TRUE
)
con <- TENxMTX(mtxf)
con
#> TENxMTX object
#> resource: /usr/local/lib/R/host-site-library/TENxIO/extdata/pbmc_3k_ff_bc_ex.mtx
import(con)
#> class: SummarizedExperiment 
#> dim: 171 10 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowData names(0):
#> colnames: NULL
#> colData names(0):
```

# TENxFileList

The `TENxFileList` class easily allows importing multiple files within a
`tar.gz` archive. The `untar` function can list the files compressed
within the tarball.

``` r
fl <- system.file(
    "extdata", "pbmc_granulocyte_sorted_3k_ff_bc_ex_matrix.tar.gz",
    package = "TENxIO", mustWork = TRUE
)
untar(fl, list = TRUE)
#> [1] "./pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix/filtered_feature_bc_matrix/"               
#> [2] "./pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
#> [3] "./pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix/filtered_feature_bc_matrix/features.tsv.gz"
#> [4] "./pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix/filtered_feature_bc_matrix/matrix.mtx.gz"
```

Using a similar import process across all file types, one can easily
obtain a Bioconductor representation that is ready for analysis.
`TENxFileList` can be imported to `SingleCellExperiment`.

``` r
con <- TENxFileList(fl)
import(con)
#> class: SingleCellExperiment 
#> dim: 10 10 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowData names(3): ID Symbol Type
#> colnames(10): AAACAGCCAAATATCC-1 AAACAGCCAGGAACTG-1 ... AAACCGCGTGAGGTAG-1 AAACGCGCATACCCGG-1
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: Gene Expression
#> altExpNames(0):
```

# TENxPeaks

Peak files can be handled with the `TENxPeaks` class. These files are
usually named `*peak_annotation` files with a `.tsv` extension. Peak
files are represented as `GRanges`.

``` r
pfl <- system.file(
    "extdata", "pbmc_granulocyte_sorted_3k_ex_atac_peak_annotation.tsv",
    package = "TENxIO", mustWork = TRUE
)
tenxp <- TENxPeaks(pfl)
peak_anno <- import(tenxp)
peak_anno
#> GRanges object with 10 ranges and 3 metadata columns:
#>        seqnames        ranges strand |        gene  distance   peak_type
#>           <Rle>     <IRanges>  <Rle> | <character> <numeric> <character>
#>    [1]     chr1    9768-10660      * | MIR1302-2HG    -18894      distal
#>    [2]     chr1 180582-181297      * |  AL627309.5     -6721      distal
#>    [3]     chr1 181404-181887      * |  AL627309.5     -7543      distal
#>    [4]     chr1 191175-192089      * |  AL627309.5    -17314      distal
#>    [5]     chr1 267561-268455      * |  AP006222.2       707      distal
#>    [6]     chr1 270864-271747      * |  AP006222.2      4010      distal
#>    [7]     chr1 273947-274758      * |  AP006222.2      7093      distal
#>    [8]     chr1 585751-586647      * |  AC114498.1      -982    promoter
#>    [9]     chr1 629484-630393      * |  AC114498.1     41856      distal
#>   [10]     chr1 633556-634476      * |  AC114498.1     45928      distal
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

# TENxFragments

Fragment files are quite large and we make use of the `Rsamtools`
package to import them with the `yieldSize` parameter. By default, we
use a `yieldSize` of 200.

``` r
fr <- system.file(
    "extdata", "pbmc_3k_atac_ex_fragments.tsv.gz",
    package = "TENxIO", mustWork = TRUE
)
```

Internally, we use the `TabixFile` constructor function to work with
indexed `tsv.gz` files.

**Note**. A warning is emitted whenever a `yieldSize` parameter is not
set.

``` r
tfr <- TENxFragments(fr)
#> Warning in TENxFragments(fr): Using default 'yieldSize' parameter
tfr
#> TENxFragments object
#> resource: /home/rstudio/gh/TENxIO/inst/extdata/pbmc_3k_atac_ex_fragments.tsv.gz
```

Because there may be a variable number of fragments per barcode, we use
a `RaggedExperiment` representation for this file type.

``` r
fra <- import(tfr)
fra
#> class: RaggedExperiment 
#> dim: 10 10 
#> assays(2): barcode readSupport
#> rownames: NULL
#> colnames(10): AAACCGCGTGAGGTAG-1 AAGCCTCCACACTAAT-1 ... TGATTAGTCTACCTGC-1 TTTAGCAAGGTAGCTT-1
#> colData names(0):
```

Similar operations to those used with SummarizedExperiment are
supported. For example, the genomic ranges can be displayed via
`rowRanges`:

``` r
rowRanges(fra)
#> GRanges object with 10 ranges and 0 metadata columns:
#>        seqnames      ranges strand
#>           <Rle>   <IRanges>  <Rle>
#>    [1]     chr1 10152-10180      *
#>    [2]     chr1 10152-10195      *
#>    [3]     chr1 10080-10333      *
#>    [4]     chr1 10091-10346      *
#>    [5]     chr1 10152-10180      *
#>    [6]     chr1 10152-10202      *
#>    [7]     chr1 10097-10344      *
#>    [8]     chr1 10080-10285      *
#>    [9]     chr1 10090-10560      *
#>   [10]     chr1 10074-10209      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

# Session Information

``` r
sessionInfo()
#> R version 4.2.1 (2022-06-23)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 20.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#>  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] TENxIO_0.99.4               ExperimentHub_2.5.0         AnnotationHub_3.5.2         BiocFileCache_2.5.0         dbplyr_2.2.1               
#>  [6] rhdf5_2.41.1                SingleCellExperiment_1.19.1 SummarizedExperiment_1.27.3 Biobase_2.57.1              GenomicRanges_1.49.1       
#> [11] GenomeInfoDb_1.33.7         IRanges_2.31.2              S4Vectors_0.35.4            BiocGenerics_0.43.4         MatrixGenerics_1.9.1       
#> [16] matrixStats_0.62.0         
#> 
#> loaded via a namespace (and not attached):
#>   [1] ellipsis_0.3.2                rprojroot_2.0.3               XVector_0.37.1                fs_1.5.2                     
#>   [5] rstudioapi_0.14               remotes_2.4.2                 bit64_4.0.5                   interactiveDisplayBase_1.35.0
#>   [9] AnnotationDbi_1.59.1          fansi_1.0.3                   codetools_0.2-18              R.methodsS3_1.8.2            
#>  [13] cachem_1.0.6                  knitr_1.40                    pkgload_1.3.0                 Rsamtools_2.13.4             
#>  [17] png_0.1-7                     R.oo_1.25.0                   shiny_1.7.2                   HDF5Array_1.25.2             
#>  [21] BiocManager_1.30.18           readr_2.1.3                   compiler_4.2.1                httr_1.4.4                   
#>  [25] assertthat_0.2.1              Matrix_1.5-1                  fastmap_1.1.0                 cli_3.4.1                    
#>  [29] later_1.3.0                   htmltools_0.5.3               prettyunits_1.1.1             tools_4.2.1                  
#>  [33] glue_1.6.2                    GenomeInfoDbData_1.2.9        dplyr_1.0.10                  rappdirs_0.3.3               
#>  [37] Rcpp_1.0.9                    vctrs_0.4.2                   Biostrings_2.65.6             rhdf5filters_1.9.0           
#>  [41] RaggedExperiment_1.21.7       xfun_0.33                     stringr_1.4.1                 ps_1.7.1                     
#>  [45] brio_1.1.3                    testthat_3.1.4                mime_0.12                     miniUI_0.1.1.1               
#>  [49] lifecycle_1.0.2               devtools_2.4.4                zlibbioc_1.43.0               vroom_1.6.0                  
#>  [53] hms_1.1.2                     promises_1.2.0.1              parallel_4.2.1                yaml_2.3.5                   
#>  [57] curl_4.3.2                    memoise_2.0.1                 stringi_1.7.8                 RSQLite_2.2.18               
#>  [61] BiocVersion_3.16.0            BiocIO_1.7.1                  desc_1.4.2                    filelock_1.0.2               
#>  [65] BiocParallel_1.31.12          pkgbuild_1.3.1                rlang_1.0.6                   pkgconfig_2.0.3              
#>  [69] bitops_1.0-7                  evaluate_0.16                 lattice_0.20-45               purrr_0.3.4                  
#>  [73] Rhdf5lib_1.19.2               htmlwidgets_1.5.4             bit_4.0.4                     processx_3.7.0               
#>  [77] tidyselect_1.1.2              magrittr_2.0.3                R6_2.5.1                      generics_0.1.3               
#>  [81] profvis_0.3.7                 DelayedArray_0.23.2           DBI_1.1.3                     pillar_1.8.1                 
#>  [85] withr_2.5.0                   colorout_1.2-2                KEGGREST_1.37.3               RCurl_1.98-1.9               
#>  [89] tibble_3.1.8                  crayon_1.5.2                  utf8_1.2.2                    tzdb_0.3.0                   
#>  [93] rmarkdown_2.16                urlchecker_1.0.1              usethis_2.1.6                 grid_4.2.1                   
#>  [97] blob_1.2.3                    callr_3.7.2                   digest_0.6.29                 xtable_1.8-4                 
#> [101] httpuv_1.6.6                  BiocBaseUtils_0.99.12         R.utils_2.12.0                sessioninfo_1.2.2
```
