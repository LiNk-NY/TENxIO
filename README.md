
<!-- badges: start -->
<!-- badges: end -->

# Introduction

`TENxIO` allows users to import 10X pipeline files into known
Bioconductor classes. The package is not comprehensive, there are files
that are not supported. It currently does not support Visium datasets.
It does replace some functionality in `DropletUtils`. If you would like
a file format to be supported. Please open an issue at
<https://github.com/waldronlab/TENxIO>.

# Supported Formats

| **Extension**       | **Class**     | **Imported as**                    |
|---------------------|---------------|------------------------------------|
| .h5                 | TENxH5        | SingleCellExperiment w/ TENxMatrix |
| .mtx / .mtx.gz      | TENxMTX       | SummarizedExperiment w/ dgCMatrix  |
| .tar.gz             | TENxFileList  | SingleCellExperiment w/ dgCMatrix  |
| peak_annotation.tsv | TENxPeaks     | GRanges                            |
| fragments.tsv.gz    | TENxFragments | RaggedExperiment                   |
| .tsv / .tsv.gz      | TSVFile\*     | tibble                             |

**Note** (\*). The `TSVFile` class is used internally and not exported.

# Tested 10X Products

We have tested these functions with *some*
[datasets](https://www.10xgenomics.com/resources/datasets) from 10x
Genomics including those from:

- Single Cell Gene Expression
- Single Cell ATAC
- Single Cell Multiome ATAC + Gene Expression

Note. That extensive testing has not been performed and the codebase may
require some adaptation to ensure compatibility with all pipeline
outputs.

# Currently not supported

- Spatial Gene Expression

# Bioconductor implementations

We are aware of existing functionality in both `DropletUtils` and
`SpatialExperiment`. We are working with the authors of those packages
to cover the use cases in both those packages and possibly port I/O
functionality into `TENxIO`. We are using long tests and the
`DropletTestFiles` package to cover example datasets on `ExperimentHub`,
if you would like to know more, see the `longtests` directory on GitHub.

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
files typically obtained from the 10X Genomics website. Generally, these
are outputs from the Cell Ranger pipeline.

# Procedure

Loading the data into a Bioconductor class is a two step process. First,
the file must be identified by either the user or the `TENxFile`
function. The appropriate function will be evoked to provide a `TENxIO`
class representation, e.g., `TENxH5` for HDF5 files with an `.h5`
extension. Secondly, the `import` method for that particular file class
will render a common Bioconductor class representation for the user. The
main representations used by the package are `SingleCellExperiment`,
`SummarizedExperiment`, `GRanges`, and `RaggedExperiment`.

# Dataset versioning

The versioning schema in the package mostly applies to HDF5 resources
and is loosely based on versions of 10X datasets. For the most part,
version 3 datasets usually contain ranged information at specific
locations in the data file. Version 2 datasets will usually contain a
`genes.tsv` file, rather than `features.tsv` as in version 3. If the
file version is unknown, the software will attempt to derive the version
from the data where possible.

# File classes

## TENxFile

The `TENxFile` class is the catch-all class superclass that allows
transition to subclasses pertinent to specific files. It inherits from
the `BiocFile` class and allows for easy dispatching `import` methods.

``` r
showClass("TENxFile")
#> Class "TENxFile" [package "TENxIO"]
#> 
#> Slots:
#>                                                                                                                                                       
#> Name:                extension                  colidx                  rowidx                  remote              compressed                resource
#> Class:               character                 integer                 integer                 logical                 logical character_OR_connection
#> 
#> Extends: "BiocFile"
#> 
#> Known Subclasses: "TSVFile", "TENxFragments", "TENxH5", "TENxMTX", "TENxPeaks"
```

### `ExperimentHub` resources

`TENxFile` can handle resources from `ExperimentHub` with careful
inputs. For example, one can import a `TENxBrainData` dataset via the
appropriate `ExperimentHub` identifier (`EH1039`):

``` r
hub <- ExperimentHub::ExperimentHub()
#> snapshotDate(): 2022-10-24
hub["EH1039"]
#> ExperimentHub with 1 record
#> # snapshotDate(): 2022-10-24
#> # names(): EH1039
#> # package(): TENxBrainData
#> # $dataprovider: 10X Genomics
#> # $species: Mus musculus
#> # $rdataclass: character
#> # $rdatadateadded: 2017-10-26
#> # $title: Brain scRNA-seq data, 'HDF5-based 10X Genomics' format
#> # $description: Single-cell RNA-seq data for 1.3 million brain cells from E18 mice. 'HDF5-based 10X Genomics' format originally provided by TENx Ge...
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

Note. `EH1039` is a large \~ 4GB file and files without extension as
those obtained from `ExperimentHub` will emit a warning so that the user
is aware that the import operation may fail, esp. if the internal
structure of the file is modified.

### TENxH5

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

Note. The `h5ls` function gives an overview of the structure of the
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

## TENxMTX

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
```

## import MTX method

The `import` method yields a `SummarizedExperiment` without colnames or
rownames.

``` r
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

## TENxFileList

Generally, the 10X website will provide tarballs (with a `.tar.gz`
extension) which can be imported with the `TENxFileList` class. The
tarball can contain components of a gene expression experiment including
the matrix data, row data (aka ‘features’) expressed as Ensembl
identifiers, gene symbols, etc. and barcode information for the columns.

The `TENxFileList` class allows importing multiple files within a
`tar.gz` archive. The `untar` function with the `list = TRUE` argument
shows all the file names in the tarball.

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

We then use the `import` method across all file types to obtain an
integrated Bioconductor representation that is ready for analysis. Files
in `TENxFileList` can be represented as a `SingleCellExperiment` with
row names and column names.

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

## TENxPeaks

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

## TENxFragments

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
#> resource: /usr/local/lib/R/host-site-library/TENxIO/extdata/pbmc_3k_atac_ex_fragments.tsv.gz
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
#>  [1] rhdf5_2.41.1                TENxIO_0.99.9               SingleCellExperiment_1.19.1 SummarizedExperiment_1.27.3 Biobase_2.57.1             
#>  [6] GenomicRanges_1.49.1        GenomeInfoDb_1.33.15        IRanges_2.31.2              S4Vectors_0.35.4            BiocGenerics_0.43.4        
#> [11] MatrixGenerics_1.9.1        matrixStats_0.62.0          colorout_1.2-2             
#> 
#> loaded via a namespace (and not attached):
#>  [1] bitops_1.0-7                  bit64_4.0.5                   filelock_1.0.2                httr_1.4.4                   
#>  [5] tools_4.2.1                   utf8_1.2.2                    R6_2.5.1                      HDF5Array_1.25.2             
#>  [9] DBI_1.1.3                     rhdf5filters_1.9.0            withr_2.5.0                   tidyselect_1.2.0             
#> [13] bit_4.0.4                     curl_4.3.3                    compiler_4.2.1                cli_3.4.1                    
#> [17] DelayedArray_0.23.2           readr_2.1.3                   rappdirs_0.3.3                Rsamtools_2.13.4             
#> [21] stringr_1.4.1                 digest_0.6.30                 rmarkdown_2.17                R.utils_2.12.0               
#> [25] XVector_0.37.1                pkgconfig_2.0.3               htmltools_0.5.3               dbplyr_2.2.1                 
#> [29] fastmap_1.1.0                 rlang_1.0.6                   rstudioapi_0.14               RSQLite_2.2.18               
#> [33] shiny_1.7.2                   BiocIO_1.7.1                  generics_0.1.3                BiocParallel_1.31.14         
#> [37] vroom_1.6.0                   dplyr_1.0.10                  R.oo_1.25.0                   RCurl_1.98-1.9               
#> [41] magrittr_2.0.3                GenomeInfoDbData_1.2.9        Matrix_1.5-1                  Rcpp_1.0.9                   
#> [45] Rhdf5lib_1.19.2               fansi_1.0.3                   lifecycle_1.0.3               R.methodsS3_1.8.2            
#> [49] stringi_1.7.8                 yaml_2.3.6                    RaggedExperiment_1.21.8       zlibbioc_1.43.0              
#> [53] BiocFileCache_2.5.2           AnnotationHub_3.5.2           grid_4.2.1                    blob_1.2.3                   
#> [57] parallel_4.2.1                promises_1.2.0.1              ExperimentHub_2.5.0           crayon_1.5.2                 
#> [61] lattice_0.20-45               Biostrings_2.65.6             hms_1.1.2                     KEGGREST_1.37.3              
#> [65] knitr_1.40                    pillar_1.8.1                  codetools_0.2-18              glue_1.6.2                   
#> [69] BiocVersion_3.16.0            evaluate_0.17                 BiocManager_1.30.18           png_0.1-7                    
#> [73] vctrs_0.4.2                   tzdb_0.3.0                    httpuv_1.6.6                  purrr_0.3.5                  
#> [77] assertthat_0.2.1              cachem_1.0.6                  xfun_0.34                     BiocBaseUtils_0.99.12        
#> [81] mime_0.12                     xtable_1.8-4                  later_1.3.0                   tibble_3.1.8                 
#> [85] AnnotationDbi_1.59.1          memoise_2.0.1                 ellipsis_0.3.2                interactiveDisplayBase_1.35.1
#> [89] BiocStyle_2.25.0
```
