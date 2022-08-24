
<!-- badges: start -->
<!-- badges: end -->

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
#> Name:                extension                  colidx                  rowidx
#> Class:               character                 integer                 integer
#>                                                       
#> Name:                   remote                resource
#> Class:                 logical character_OR_connection
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
#> snapshotDate(): 2022-08-23
hub["EH1039"]
#> ExperimentHub with 1 record
#> # snapshotDate(): 2022-08-23
#> # names(): EH1039
#> # package(): TENxBrainData
#> # $dataprovider: 10X Genomics
#> # $species: Mus musculus
#> # $rdataclass: character
#> # $rdatadateadded: 2017-10-26
#> # $title: Brain scRNA-seq data, 'HDF5-based 10X Genomics' format
#> # $description: Single-cell RNA-seq data for 1.3 million brain cells from E1...
#> # $taxonomyid: 10090
#> # $genome: mm10
#> # $sourcetype: HDF5
#> # $sourceurl: http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M...
#> # $sourcesize: NA
#> # $tags: c("SequencingData", "RNASeqData", "ExpressionData",
#> #   "SingleCell") 
#> # retrieve record with 'object[["EH1039"]]'
```

Currently, `ExperimentHub` resources do not have an extension and it is
best to provide that to the `TENxFile` constructor function.

``` r
fname <- hub[["EH1039"]]
#> see ?TENxBrainData and browseVignettes('TENxBrainData') for documentation
#> loading from cache
TENxFile(fname, extension = "h5", group = "mm10", version = "2")
#> Warning: 'group' not in known 10X groups: matrix, outs
#> TENxH5 object 
#> resource: /home/mr148/.cache/R/ExperimentHub/24076533f113_1039 
#> projection: SingleCellExperiment 
#> dim: 27998 1306127 
#> rownames: ENSMUSG00000051951 ENSMUSG00000089699 ... ENSMUSG00000096730 ENSMUSG00000095742 
#> rowData names(3): ID Symbol Type 
#>   Type: ENSMUSG00000025900 ENSMUSG00000025902 ... ENSMUSG00000102343 ENSMUSG00000109048 
#> colnames: AAACCTGAGATAGGAG-1 AAACCTGAGCGGCTTC-1 ... TTTGTCAGTTAAAGTG-133 TTTGTCATCTGAAAGA-133
```

## TENxH5

One of the main data formats provided by the 10X website are HDF5 files.
To import those files, we use the TENxH5 contructor function.

``` r
h5f <- system.file(
    "extdata", "pbmc_granulocyte_ff_bc_ex.h5",
    package = "TENxIO", mustWork = TRUE
) 
h5f
#> [1] "/media/mr148/1D24A0EA4286043C/bioc-devel/TENxIO/extdata/pbmc_granulocyte_ff_bc_ex.h5"
```

Here we provide a bespoke show method for such files given that we have
some idea of the structure to expect within these files:

``` r
con <- TENxH5(h5f)
con
#> TENxH5 object 
#> resource: /media/mr148/1D24A0EA4286043C/bioc-devel/TENxIO/extdata/pbmc_granulocyte_ff_bc_ex.h5 
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
#> rownames(10): ENSG00000243485 ENSG00000237613 ... ENSG00000286448
#>   ENSG00000236601
#> rowData names(3): ID Symbol Type
#> colnames(10): AAACAGCCAAATATCC-1 AAACAGCCAGGAACTG-1 ...
#>   AAACCGCGTGAGGTAG-1 AAACGCGCATACCCGG-1
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: Gene Expression
#> altExpNames(0):
```

**Note**. Future versions of the package could support alternative
representations to `SingleCellExperiment`

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
#> resource: /media/mr148/1D24A0EA4286043C/bioc-devel/TENxIO/extdata/pbmc_3k_ff_bc_ex.mtx
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

Using a similar import process accross all file types, one can easily
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
#> colnames(10): AAACAGCCAAATATCC-1 AAACAGCCAGGAACTG-1 ...
#>   AAACCGCGTGAGGTAG-1 AAACGCGCATACCCGG-1
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
#> resource: /media/mr148/1D24A0EA4286043C/bioc-devel/TENxIO/extdata/pbmc_3k_atac_ex_fragments.tsv.gz
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
#> colnames(10): AAACCGCGTGAGGTAG-1 AAGCCTCCACACTAAT-1 ...
#>   TGATTAGTCTACCTGC-1 TTTAGCAAGGTAGCTT-1
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
#> R version 4.2.1 Patched (2022-07-29 r82650)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 20.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] TENxBrainData_1.17.0        HDF5Array_1.25.2           
#>  [3] rhdf5_2.41.1                DelayedArray_0.23.1        
#>  [5] Matrix_1.4-1                TENxIO_0.13.0              
#>  [7] SingleCellExperiment_1.19.0 SummarizedExperiment_1.27.2
#>  [9] Biobase_2.57.1              GenomicRanges_1.49.1       
#> [11] GenomeInfoDb_1.33.5         IRanges_2.31.2             
#> [13] S4Vectors_0.35.1            BiocGenerics_0.43.1        
#> [15] MatrixGenerics_1.9.1        matrixStats_0.62.0         
#> 
#> loaded via a namespace (and not attached):
#>  [1] bitops_1.0-7                  bit64_4.0.5                  
#>  [3] filelock_1.0.2                httr_1.4.4                   
#>  [5] tools_4.2.1                   utf8_1.2.2                   
#>  [7] R6_2.5.1                      DBI_1.1.3                    
#>  [9] rhdf5filters_1.9.0            tidyselect_1.1.2             
#> [11] bit_4.0.4                     curl_4.3.2                   
#> [13] compiler_4.2.1                cli_3.3.0                    
#> [15] readr_2.1.2                   rappdirs_0.3.3               
#> [17] Rsamtools_2.13.3              stringr_1.4.0                
#> [19] digest_0.6.29                 rmarkdown_2.15               
#> [21] R.utils_2.12.0                XVector_0.37.0               
#> [23] pkgconfig_2.0.3               htmltools_0.5.3              
#> [25] dbplyr_2.2.1                  fastmap_1.1.0                
#> [27] rlang_1.0.4                   rstudioapi_0.13              
#> [29] RSQLite_2.2.16                shiny_1.7.2                  
#> [31] BiocIO_1.7.1                  generics_0.1.3               
#> [33] BiocParallel_1.31.12          vroom_1.5.7                  
#> [35] dplyr_1.0.9                   R.oo_1.25.0                  
#> [37] RCurl_1.98-1.8                magrittr_2.0.3               
#> [39] GenomeInfoDbData_1.2.8        Rcpp_1.0.9                   
#> [41] Rhdf5lib_1.19.2               fansi_1.0.3                  
#> [43] lifecycle_1.0.1               R.methodsS3_1.8.2            
#> [45] stringi_1.7.8                 yaml_2.3.5                   
#> [47] RaggedExperiment_1.21.1       zlibbioc_1.43.0              
#> [49] BiocFileCache_2.5.0           AnnotationHub_3.5.0          
#> [51] grid_4.2.1                    blob_1.2.3                   
#> [53] parallel_4.2.1                promises_1.2.0.1             
#> [55] ExperimentHub_2.5.0           crayon_1.5.1                 
#> [57] lattice_0.20-45               Biostrings_2.65.2            
#> [59] hms_1.1.1                     KEGGREST_1.37.3              
#> [61] knitr_1.39                    pillar_1.8.0                 
#> [63] codetools_0.2-18              glue_1.6.2                   
#> [65] BiocVersion_3.16.0            evaluate_0.16                
#> [67] BiocManager_1.30.18           png_0.1-7                    
#> [69] vctrs_0.4.1                   tzdb_0.3.0                   
#> [71] httpuv_1.6.5                  purrr_0.3.4                  
#> [73] assertthat_0.2.1              cachem_1.0.6                 
#> [75] xfun_0.32                     BiocBaseUtils_0.99.10        
#> [77] mime_0.12                     xtable_1.8-4                 
#> [79] later_1.3.0                   tibble_3.1.8                 
#> [81] AnnotationDbi_1.59.1          memoise_2.0.1                
#> [83] ellipsis_0.3.2                interactiveDisplayBase_1.35.0
```
