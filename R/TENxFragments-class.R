#' @exportClass TENxFragments
.TENxFragments <- setClass(
    Class = "TENxFragments",
    contains = "TENxFile",
    slots = c(
        dataType = "character", roi = "GenomicRanges", yieldSize = "numeric"
    )
)

.check_fragments <- function(object) {
    fpath <- path(object)
    if (!endsWith(fpath, "_fragments.tsv.gz"))
        "Provide a 10X fragments file ending in '_fragments.tsv.gz'"
    else
        TRUE
}

.validTENxFragments <- function(object) {
    .check_fragments(object)
}

S4Vectors::setValidity2("TENxFragments", .validTENxFragments)

#' @rdname TENxFragments
#' @title Import fragments files from 10X
#' @aliases TENxFragments-class
#'
#' @param resource character(1) The file path to the fragments resource, usually
#'   a compressed tabix file with extension `.tsv.gz`.
#'
#' @param roi GRanges() A GRanges indicating the regions of interest. This
#'   get sent to `RSamtools` as the `param` input.
#'
#' @param yieldSize numeric() The number of records to read by default, 200
#'   records will be imported. A warning will be emmitted if not modified.
#'
#' @param ... Further arguments to the class generator function (currently not
#'   used)
#'
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom S4Vectors splitAsList mcols
#'
#' @examples
#'
#' con <- TENxFile("~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5")
#' con
#' sce <- import(con)
#'
#' fr <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
#' tfr <- TENxFragments(fr)
#' fra <- import(tfr)
#'
#' @export
TENxFragments <- function(resource, yieldSize = 200, roi = GRanges(), ...) {
    if (missing(yieldSize))
        warning("Using default 'yieldSize' parameter")
    if (!is(roi, "GRanges"))
        stop("'roi' input must be 'GenomicRanges'")
    type <- tail(strsplit(resource, "_|_fragments.tsv.gz")[[1]], 1L)
    .TENxFragments(
        resource = resource, dataType = type, yieldSize = yieldSize, roi = roi,
        ...
    )
}

#' @export
setMethod("import", "TENxFragments", function(con, format, text, ...) {
    .checkPkgsAvail(c("Rsamtools", "RaggedExperiment"))
    roi <- con@roi
    yieldSize <- con@yieldSize
    tb <- Rsamtools::TabixFile(path(con), yieldSize = yieldSize)
    ex <- read.table(
        textConnection(Rsamtools::scanTabix(tb, param = roi)[[1]])
    )
    ## https://support.10xgenomics.com/single-cell-atac/
    ## software/pipelines/latest/output/fragments
    names(ex) <- c("chrom", "chromStart", "chromEnd", "barcode", "readSupport")
    ggr <- makeGRangesFromDataFrame(
        ex, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
    )
    grs <- splitAsList(ggr, mcols(ggr)$barcode)
    ggrl <- as(grs, "GRangesList")
    RaggedExperiment::RaggedExperiment(ggrl)
})
