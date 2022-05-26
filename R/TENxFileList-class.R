#' @exportClass TENxFileList
.TENxFileList <- setClass(
    Class = "TENxFileList",
    contains = "SimpleList",
    slots = c(compressed = "logical")
)

.validTENxFileList <- function(object) {
    if (!object@compressed) {
        validClasses <- vapply(object, function(ifile) {
            is(ifile, "TENxFile")
        }, logical(1L))
        if (all(validClasses))
            TRUE
        else
            "Some files are not of class 'TENxFile'"
    }
}

S4Vectors::setValidity2("TENxFileList", .validTENxFileList)

#' @export
TENxFileList <- function(..., compressed = FALSE) {
    dots <- list(...)
    if (identical(length(dots[[1]]), 1L) && is.character(dots[[1]]))
        exts <- .get_ext(dots[[1]])
    else if (is.list(dots[[1]]))
        exts <- vapply(dots[[1]], .get_ext, character(1L))
    ## TODO: handle folder and list of files, in addition to tar.gz
    if (identical(exts, "tar.gz"))
        compressed <- TRUE
    asl <- SimpleList(dots[[1]])
    .TENxFileList(asl, compressed = compressed)
}
