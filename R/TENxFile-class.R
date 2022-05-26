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

#' @exportClass TENxH5
.TENxH5 <- setClass(
    Class = "TENxH5",
    contains = "TENxFile",
    slots = c(version = "character", group = "character")
)

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

.get_h5_group <- function(fpath) {
    if (!requireNamespace("rhdf5", quietly = TRUE))
        stop("Install 'rhdf5' to work with TENxH5")
    l1 <- rhdf5::h5ls(fpath, recursive = FALSE)
    l1[l1$otype == "H5I_GROUP", "name"]
}

.KNOWN_H5_GROUPS <- c("matrix", "outs")

.check_group_ok <- function(fpath) {
    gname <- .get_h5_group(fpath)
    if (!gname %in% .KNOWN_H5_GROUPS)
        stop("'group' not recognized")
    gname
}

#' @rdname TENxH5
#' @title Import H5 files from 10X
#'
#' @examples
#'
#' h5f <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
#' con <- TENxFile(h5f)
#' import(con)
#'
#' @export
TENxH5 <-
    function(resource, version = c("3", "2"), ...)
{
    version <- match.arg(version)
    group <- .check_group_ok(resource)
    .TENxH5(resource = resource, group = group, version = version, ...)
}
