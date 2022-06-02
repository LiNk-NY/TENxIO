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
        vlist <- list(.proc_peak_value(value))
        names(vlist) <- args[["name"]]
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
#' @examples
#'
#' con <- TENxFile("~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5")
#' sce <- import(con)
#' fr <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
#' fra <- .import_fragments(fr)
#'
.import_fragments <- function(file, withIndex = TRUE) {
    ## read tabix format with Rsamtools
    stopifnot(endsWith(file, "atac_fragments.tsv.gz"))
    tabindex <- paste0(file, ".tbi")
    if (withIndex)
        stopifnot(file.exists(tabindex))

    tb <- Rsamtools::TabixFile(file, tabindex, yieldSize = NA_integer_)
    # reslist <- scanTabix(tb, param = rowRanges(sce))
    # lapply(Filter(length, reslist), function(x) read.table(textConnection(x)))
    ex <- read.table(textConnection(scanTabix(tb)[[1]]))

    ## names from table in the link below
    ## https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
    names(ex) <- c("chrom", "chromStart", "chromEnd", "barcode", "readSupport")
    ggr <- makeGRangesFromDataFrame(
        ex, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
    )
    ggrl <- as(splitAsList(ggr, mcols(ggr)$barcode), "GRangesList")
    RaggedExperiment(ggrl)
    # compactAssay(erg, i = "readSupport", sparse = TRUE)
}
