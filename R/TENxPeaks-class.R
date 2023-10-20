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
#' @return A `TENxPeaks` class object
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
#' @inheritParams TENxFile
#'
#' @return A `GRanges` class object of peak locations
#'
#' @examples
#'
#' fi <- system.file(
#'     "extdata", "pbmc_granulocyte_sorted_3k_ex_atac_peak_annotation.tsv",
#'     package = "TENxIO", mustWork = TRUE
#' )
#' peak_file <- TENxPeaks(fi)
#' peak_anno <- import(peak_file)
#' peak_anno
#'
#' example(TENxH5)
#'
#' ## Add peaks to an existing SCE
#' ## First, import the SCE from an example H5 file
#' h5f <- system.file(
#'     "extdata", "pbmc_granulocyte_ff_bc_ex.h5",
#'     package = "TENxIO", mustWork = TRUE
#' )
#' con <- TENxH5(h5f)
#' sce <- import(con)
#' ## auto-import peaks when using annotation<-
#' annotation(sce, name = "peak_annotation") <- peak_file
#' annotation(sce)
#'
#' @export
TENxPeaks <- function(resource, extension, ...) {
    if (missing(extension))
        extension <- .get_ext(resource)
    .TENxPeaks(resource = resource, extension = extension, ...)
}

#' @describeIn TENxPeaks-class Import a peaks_annotation file from 10x as a
#'   `GRanges` representation
#'
#' @inheritParams BiocIO::import
#'
#' @export
setMethod("import", "TENxPeaks", function(con, format, ...) {
    .checkPkgsAvail("readr")
    if (identical(Sys.info()[["sysname"]], "Darwin"))
        readr::local_edition(1)
    panno <- readr::read_tsv(
        file = path(con), col_types = c("c", "n", "n", "c", "n", "c")
    )
    makeGRangesFromDataFrame(panno, keep.extra.columns = TRUE)
})

#' @describeIn TENxPeaks-class Replacement method to add annotation data to a
#'   `SingleCellExperiment`
#'
#' @inheritParams BiocGenerics::annotation
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom BiocGenerics annotation<-
#'
#' @export
setReplaceMethod("annotation", "SingleCellExperiment",
    function(object, ..., value) {
        if (!is(value, "TENxPeaks"))
            stop("'value' must be of class 'TENxPeaks'")
        args <- list(...)
        append <- isTRUE(args[["append"]])
        anno_name <- args[["name"]]
        if (is.null(anno_name))
            anno_name <- "peak_annotation"
        meta <- metadata(object)
        annotated <- "annotation" %in% names(meta)
        if (annotated && !append)
            warning(
                "'append = FALSE'; replacing annotation in metadata",
                call. = FALSE
            )
        vlist <- structure(
            list(import(value)), .Names = anno_name
        )
        if (append)
            vlist <- append(meta[["annotation"]], vlist)
        metadata(object)[["annotation"]] <- vlist
        object
    }
)

#' @describeIn TENxPeaks-class Extraction method to obtain annotation data from
#'   a `SingleCellExperiment` representation
#'
#' @importFrom BiocGenerics annotation
#'
#' @export
setMethod("annotation", "SingleCellExperiment", function(object, ...) {
    metadata(object)[["annotation"]]
})
