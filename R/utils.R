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
    .checkPkgsAvail("readr")
    panno <- readr::read_tsv(
        file = file, col_types = c("c", "n", "n", "c", "n", "c")
    )
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

.checkPkgsAvail <- function(pkgnames, walkback = -3L) {
    vapply(pkgnames, function(pkgname) {
        func <- as.character(sys.call(walkback)[[1L]])
        func <- tail(func, 1L)
        if (!requireNamespace(pkgname, quietly = TRUE))
            stop("Install '", pkgname, "' to use '", func, "'", call. = FALSE)
        else
            TRUE
    }, logical(1L))
}
