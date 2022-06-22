#' TENxMTX: The Matrix Market representation class for 10X Data
#'
#'
.TENxMTX <- setClass(
    Class = "TENxMTX",
    contains = "TENxFile",
    slots = c(compressed = "logical")
)

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

#' @import SummarizedExperiment
#' @export
setMethod("import", "TENxMTX", function(con, format, text, ...) {
    mtxf <- Matrix::readMM(path(con))
    ## coerce to common use class
    mtxf <- as(mtxf, "dgCMatrix")
    SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(counts = mtxf)
    )
})
