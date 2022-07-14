#' TENxPeaks: The class to represent 10x Peaks files
#'
#' This class is designed to work with the files denoted with "peak_annotation"
#' in the file name. These are usually produced as tab separated value files,
#' i.e., `.tsv`.
#'
#' @details This class is a straightforward class for handling peak data. It can
#'   be used in conjunction with the `annotation` method on a
#'   `SingleCellExperiment` to add peak information to the experiment. The
#'   ranged data is represented as a `GRanges` class object.
#'
#' @exportClass TENxPeaks
.TENxPeaks <- setClass(
    Class = "TENxPeaks",
    contains = "TENxFile"
)

.check_peaks <- function(object) {
    fpath <- path(object)
    if (!endsWith(fpath, "peak_annotation.tsv"))
        "Provide a 10X Peaks file ending in 'peak_annotation.tsv'"
    else
        TRUE
}

.validPeaksFile <- function(object) {
    .check_peaks(object)
}

S4Vectors::setValidity2("TENxPeaks", .validPeaksFile)

#' Import 10x peak annotation files from 10x
#'
#' This constructor function is designed to work with the files denoted with
#' "peak_annotation" in the file name. These are usually produced as tab
#' separated value files, i.e., `.tsv`.
#'
#' @details The output class allows handling of peak data. It can be used in
#'   conjunction with the `annotation` method on a `SingleCellExperiment` to add
#'   peak information to the experiment. The ranged data is represented as a
#'   `GRanges` class object.
#'
#' @return A `TENxPeak` class object
#'
#' @examples
#'
#' fi <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_atac_peak_annotation.tsv"
#' peak_file <- TENxPeaks(fi)
#' peak_anno <- import(pa)
#'
#' con <- TENxFile("~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5")
#' sce <- import(con)
#' annotation(sce, name = "peak_annotation") <- peak_file
#'
TENxPeaks <- function(resource, ...) {
    .TENxPeaks(resource = resource, ...)
}

#' @export
setMethod("import", "TENxPeaks", function(con, format, ...) {
    .checkPkgsAvail("readr")
    panno <- readr::read_tsv(
        file = path(con), col_types = c("c", "n", "n", "c", "n", "c")
    )
    makeGRangesFromDataFrame(panno, keep.extra.columns = TRUE)
})

#' @importFrom S4Vectors metadata metadata<-
setReplaceMethod("annotation", "SingleCellExperiment",
    function(object, ..., value) {
        if (!is(value, "TENxPeaks"))
            stop("'value' must be of class 'TENxPeaks'")
        args <- list(...)
        append <- isTRUE(args[["append"]])
        meta <- metadata(object)
        annotated <- "annotation" %in% names(meta)
        if (annotated && !append)
            warning("'append = FALSE'; replacing annotation in metadata")
        vlist <- structure(
            list(import(value)), .Names = args[["name"]]
        )
        metadata(object) <- append(meta, list(annotation = vlist))
        object
    }
)

#' @export
setMethod("annotation", "SingleCellExperiment", function(object, ...) {
    metadata(object)[["annotation"]]
})
