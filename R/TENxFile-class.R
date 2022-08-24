#' TENxFile: General purpose class for 10X files
#'
#' @description The `TENxFile` class is the default representation for
#'   unrecognized subclasses. It inherits from the BiocFile class and adds a few
#'   additional slots. The constructor function can handle typical 10X file
#'   types. For more details, see the constructor function documentation.
#'
#' @slot extension character(1) The file extension as extracted from the file
#'   path or overridden via the `ext` argument in the constructor function.
#'
#' @slot colidx integer(1) The column index corresponding to the columns in the
#'   file that will subsequently be imported
#'
#' @slot rowidx integer(1) The row index corresponding to rows in the file that
#'   will subsequently be imported
#'
#' @slot remote logical(1) Whether the file exists on the web, i.e., the
#'   `resource` is a URL
#'
#' @importClassesFrom BiocIO BiocFile
#' @importFrom BiocIO import
#' @importFrom BiocGenerics path
#' @importFrom S4Vectors mcols
#' @importFrom methods new is as
#'
#' @exportClass TENxFile
.TENxFile <- setClass(
    Class = "TENxFile",
    contains = "BiocFile",
    slots = c(
        extension = "character",
        colidx = "integer", rowidx = "integer",
        remote = "logical"
    )
)

.check_file_exists <- function(object) {
    op <- path(object)
    if (file.exists(op) || RCurl::url.exists(op))
        TRUE
    else
        "Path or URL to the file must be valid"
}

.validTENxFile <- function(object) {
    .check_file_exists(object)
}

S4Vectors::setValidity2("TENxFile", .validTENxFile)

# TENxFile constructor ----------------------------------------------------

.remove_query <- function(fname) {
    vapply(strsplit(fname, "\\?"), `[[`, character(1L), 1L)
}

.get_ext <- function(fname) {
    if (is(fname, "TENxFile"))
        fname@extension
    else {
        fname <- .remove_query(fname)
        split_files <- strsplit(basename(fname), "\\.")
        vapply(split_files, function(file) {
            paste0(
                utils::tail(file, -1),
                collapse = "."
            )
        }, character(1L))
    }
}

#' TENxFile constructor function
#'
#' @description The `TENxFile` constructor function serves as the
#'   auto-recognizer function for 10X files. It can import several different
#'   file extensions, namely:
#' \preformatted{
#'     * H5 - on-disk HDF5
#'     * MTX - matrix market
#'     * .tar.gz - compressed tarball
#' }
#'
#' @param resource character(1) The path to the file
#'
#' @param extension character(1) The file extension for the given resource. It
#'   can usually be obtained from the file path. An override can be provided
#'   especially for `ExperimentHub` resources where the file extension is
#'   removed.
#'
#' @param ... Additional inputs to the low level class generator functions
#'
#' @return A subclass of `TENxFile` according to the input file extension
#'
#' @examples
#'
#' ## from ExperimentHub
#' hub <- ExperimentHub::ExperimentHub()
#' fname <- hub[["EH1039"]]
#' TENxFile(fname, extension = "h5", group = "mm10", version = "2")
#'
#' @export
TENxFile <- function(resource, extension, ...) {
    if (missing(extension))
        extension <- .get_ext(resource)
    TENxFUN <- switch(
        extension,
        h5 = TENxH5,
        mtx = TENxMTX,
        mtx.gz = TENxMTX,
        tar.gz = TENxFileList,
        tsv.gz = .TSVFile,
        .TENxFile
    )
    TENxFUN(resource = resource,  extension = extension, ...)
}
