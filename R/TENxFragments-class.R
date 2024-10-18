#' @include TENxFile-class.R
#'
NULL

#' TENxFragments: A class to represent fragments data as `GRanges`
#'
#' This class is designed to work mainly with `fragments.tsv.gz` files from
#' 10x pipelines.
#'
#' @details Fragments data from 10x can be quite large. In order to speed up
#' the initial exploration of the data, we use a default of **200** records
#' for loading. Users can change this default value by specifying a new one
#' via the `yieldSize` argument in the constructor function.
#'
#' @slot which GRanges() A GRanges indicating the regions of interest. This
#'   get sent to `RSamtools` as the `param` input.
#'
#' @slot yieldSize numeric() The number of records to read by default, 200
#'   records will be imported. A warning will be emitted if not modified.
#'
#' @return A `TENxFragments` class object
#'
#' @exportClass TENxFragments
.TENxFragments <- setClass(
    Class = "TENxFragments",
    contains = "TENxFile",
    slots = c(which = "GenomicRanges", yieldSize = "numeric")
)

.check_fragments <- function(object) {
    fpath <- path(object)
    fragex <- endsWith(fpath, c("_fragments.tsv.gz", "_fragments.tsv"))
    if (!any(fragex))
        "Provide a 10X fragments file ending in '_fragments.tsv*'"
    else
        TRUE
}

.validTENxFragments <- function(object) {
    .check_fragments(object)
}

S4Vectors::setValidity2("TENxFragments", .validTENxFragments)

#' TENxFragments: Import fragments files from 10X
#'
#' @param resource character(1) The file path to the fragments resource, usually
#'   a compressed tabix file with extension `.tsv.gz`.
#'
#' @param which GRanges() A GRanges indicating the regions of interest. This
#'   get sent to `RSamtools` as the `param` input.
#'
#' @param yieldSize numeric() The number of records to read by default, 200
#'   records will be imported. A warning will be emitted if not modified.
#'
#' @param ... Further arguments to the class generator function (currently not
#'   used)
#'
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom S4Vectors splitAsList mcols
#'
#' @return A `RaggedExperiment` object class
#'
#' @examples
#'
#' fr <- system.file(
#'     "extdata", "pbmc_3k_atac_ex_fragments.tsv.gz",
#'     package = "TENxIO", mustWork = TRUE
#' )
#'
#' tfr <- TENxFragments(fr)
#'
#' fra <- import(tfr)
#'
#' @export
TENxFragments <- function(resource, yieldSize = 200, which = GRanges(), ...) {
    if (missing(yieldSize) && missing(which))
        warning("Using default 'yieldSize' parameter")
    else if (!missing(which))
        yieldSize <- NA_integer_
    if (!is(which, "GRanges"))
        stop("'which' input must be 'GenomicRanges'")
    .TENxFragments(
        resource = resource, yieldSize = yieldSize, which = which,
        ...
    )
}

#' @describeIn TENxFragments Import method for representing fragments.tsv.gz
#'   data from 10x via `Rsamtools` and `RaggedExperiment`
#'
#' @importFrom utils read.table
#' @importFrom BiocBaseUtils checkInstalled
#'
#' @inheritParams BiocIO::import
#'
#' @export
setMethod("import", "TENxFragments", function(con, format, text, ...) {
    checkInstalled(c("Rsamtools", "RaggedExperiment"))
    which <- con@which
    yieldSize <- con@yieldSize
    tb <- Rsamtools::TabixFile(path(con), yieldSize = yieldSize)
    ex <- read.table(
        textConnection(Rsamtools::scanTabix(tb, param = which)[[1]])
    )
    ## https://support.10xgenomics.com/single-cell-atac/
    ## software/pipelines/latest/output/fragments
    names(ex) <- c("chrom", "chromStart", "chromEnd", "barcode", "readSupport")
    ggr <- makeGRangesFromDataFrame(
        ex, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
    )
    grs <- splitAsList(ggr, mcols(ggr)$barcode)
    ggrl <- as(grs, "GRangesList")
    RaggedExperiment::RaggedExperiment(ggrl, metadata = metadata(con))
})
