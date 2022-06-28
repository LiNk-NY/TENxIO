#' TENxIO: A Bioconductor package for importing 10X Genomics files
#'
#' The package provides file classes based on `BiocIO` for common file
#' extensions found in the 10X Genomics website.
#'
#' @section Supported file types:
#' Here is a table of supported file and file extensions and their imported
#' classes:
#'
#' | **Extension** | **Class** | **Imported as** |
#' |---------------|-----------|-----------------|
#' | .h5 | TENxH5 | SingleCellExperiment |
#' | .mtx / .mtx.gz | TENxMTX | dgCMatrix |
#' | .tar.gz | TENxFileList | SingleCellExperiment |
#' | peak_annotation.tsv | TENxPeaks | GRanges |
#' | fragments.tsv.gz | TENxFragments | RaggedExperiment |
#' | .tsv / .tsv.gz | TSVFile* | tibble |
#'
#' Note (*). The `TSVFile` class is used internally and not exported.
#'
#' @import SummarizedExperiment SingleCellExperiment
#'
#' @docType package
#'
#' @aliases TENxIO-package
#'
#' @name TENxIO
#'
"_PACKAGE"
