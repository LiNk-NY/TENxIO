#' @importClassesFrom BiocIO BiocFile
#' @importFrom BiocIO import
#' @importFrom BiocGenerics path
#' @importFrom S4Vectors mcols
#'
#' @exportClass TENxFile
.TENxFile <- setClass(
    Class = "TENxFile",
    contains = "BiocFile",
    slots = c(extension = "character")
)

.check_file_exists <- function(object) {
    if (file.exists(path(object)))
        TRUE
    else
        "Path to the file must be valid"
}

.validTENxFile <- function(object) {
    .check_file_exists(object)
}

S4Vectors::setValidity2("TENxFile", .validTENxFile)

.TENxMTX <- setClass(
    Class = "TENxMTX",
    contains = "TENxFile"
)

# TENxFile constructor ----------------------------------------------------

.get_ext <- function(fname) {
    split_files <- strsplit(basename(fname), "\\.")
    vapply(split_files, function(file) {
        paste0(
            utils::tail(file, -1),
            collapse = "."
        )
    }, character(1L))
}

#' TENxFile constructor function
#'
#' @examples
#'
#' h5f <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
#' con <- TENxFile(h5f)
#' import(con)
#'
#' ## compressed
#' h5c <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.tar.gz"
#' comp <- TENxFile(h5c)
#' import(comp)
#'
#' @export
TENxFile <- function(resource, ...) {
    ext <- .get_ext(resource)
    TENxFUN <- switch(
        ext,
        h5 = TENxH5, mtx = .TENxMTX, tar.gz = .TENxCompressed, .TENxFile
    )
    TENxFUN(resource = resource,  extension = ext, ...)
}

