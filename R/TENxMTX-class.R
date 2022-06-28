#' TENxMTX: The Matrix Market representation class for 10X Data
#'
#' @description This class is designed to work with 10x MTX datasets,
#'   particularly from the multiome pipelines.
#'
#' @slot compressed logical(1) Whether or not the file is in compressed format,
#'   usually gzipped (`.gz`).
#'
#' @details The `TENxMTX` class is a straightforward implementation that allows
#'   the user to import a Matrix Market file format using `Matrix::readMM`.
#'   Currently, it only supports return types of `dgCMatrix`. To request other
#'   formats, please open an issue on GitHub.
#'
#' @exportClass TENxMTX
.TENxMTX <- setClass(
    Class = "TENxMTX",
    contains = "TENxFile",
    slots = c(compressed = "logical")
)

#' TENxMTX: Represent Matrix Market Format Files from 10X
#'
#' This constructor function accepts `.mtx` and `.mtx.gz` compressed formats
#' for eventual importing. It is mainly used with tarball files from 10X
#' Genomics, where more annotation data is included. Importing solely the
#' `.mtx` format will provide users with a sparse matrix of `dgCMatrix` class
#' from the `Matrix` package. Currently, other formats are not supported but
#' if you'd like to request support for a format, please open an issue on
#' GitHub.
#'
#' @inheritParams TENxFile
#'
#' @param compressed logical(1) Whether the resource file is compressed (default
#'   FALSE)
#'
#' @return An instance of the `TENxMTX` class
#'
#' @examples
#'
#' mtxf <-"~/data/10x/pbmc_3k/filtered_feature_bc_matrix/matrix.mtx.gz"
#' con <- TENxMTX(mtxf)
#' import(con)
#'
#' @export
TENxMTX <- function(resource, compressed = FALSE, ...) {
    dots <- list(...)
    ext <- dots[["extension"]]
    if (is.null(ext))
        ext <- .get_ext(resource)
    compr <- identical(ext, "mtx.gz")
    if (!ext %in% c("mtx.gz", "mtx"))
        warning("File extension is not 'mtx'; import may fail", call. = FALSE)
    .TENxMTX(resource = resource, compressed = compr, extension = ext)
}

#' @export
setMethod("import", "TENxMTX", function(con, format, text, ...) {
    mtxf <- Matrix::readMM(path(con))
    ## coerce to common use class
    mtxf <- as(mtxf, "dgCMatrix")
    SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(counts = mtxf)
    )
})
