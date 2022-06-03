.setNames <- function(object = nm, nm) {
    names(object) <- nm
    object
}

#' Importing peak-annotation files
#'
#' @keywords internal
#'
#' @examples
#'
#' fi <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_atac_peak_annotation.tsv"
#' pa <- .import_peak_anno(fi)
#'
#' con <- TENxFile("~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5")
#' sce <- import(con)
#' annotation(sce, name = "peak_annotation") <- pa
#'
.import_peak_anno <- function(file) {
    stopifnot(endsWith(file, "peak_annotation.tsv"))
    panno <- readr::read_tsv(file, col_types = c("c", "n", "n", "c", "n", "c"))
    makeGRangesFromDataFrame(panno, keep.extra.columns = TRUE)
}

.proc_peak_value <- function(value) {
    if (is.character(value) && file.exists(value))
        value <- .import_peak_anno(value)
    value
}

setReplaceMethod("annotation", "SingleCellExperiment",
    function(object, ..., value) {
        args <- list(...)
        append <- isTRUE(args[["append"]])
        meta <- metadata(object)
        annotated <- "annotation" %in% names(meta)
        if (annotated && !append)
            warning("'append = FALSE'; replacing annotation in metadata")
        vlist <- structure(
            list(.proc_peak_value(value)), .Names = args[["name"]]
        )
        metadata(object) <- append(meta, list(annotation = vlist))
        object
    }
)

setMethod("annotation", "SingleCellExperiment", function(object, ...) {
    metadata(object)[["annotation"]]
})

#' Importing fragment files
#'
#' @keywords internal
#'
#' @param file character(1) The file path to the fragments resource, usually a
#'   compressed tabix file with extension `.tsv.gz`.
#'
#' @param roi GRanges() A GRanges indicating the regions of interest. This
#'   get sent to `RSamtools` as the `param` input.
#'
#' @param yieldSize numeric() The number of records to read by default, 200
#'   records will be imported. A warning will be emmitted if not modified.
#'
#' @param index character(1) The index file usually ending in `.tbi`. It
#'   should be found in the same location as the `file` input
#'
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom S4Vectors splitAsList mcols
#'
#' @examples
#'
#' con <- TENxFile("~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5")
#' sce <- import(con)
#'
#' fr <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
#' fra <- .import_fragments(fr)
#'
.import_fragments <-
    function(
        file, roi = GRanges(), yieldSize = 200, index = paste0(file, ".tbi")
    )
{
    .checkPkgsAvail(c("Rsamtools", "RaggedExperiment"))
    if (!endsWith(file, "atac_fragments.tsv.gz"))
        stop("Provide a 10X fragments file ending in 'atac_fragments.tsv.gz'")

    if (!file.exists(index))
        stop("Tabix file index not found")

    if (missing(yieldSize))
        warning("Using default 'yieldSize' parameter")
    tb <- Rsamtools::TabixFile(file, tabindex, yieldSize = yieldSize)
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
}

.checkPkgsAvail <- function(pkgnames) {
    vapply(pkgnames, function(pkgname) {
        func <- as.character(sys.call(-3L)[[1L]])
        func <- tail(func, 1L)
        if (!requireNamespace(pkgname, quietly = TRUE))
            stop("Install '", pkgname, "' to use '", func, "'", call. = FALSE)
        else
            TRUE
    }, logical(1L))
}
