#' @importClassesFrom BiocIO BiocFile
#' @importFrom BiocIO import
#' @importFrom BiocGenerics path
#' @importFrom S4Vectors mcols
#'
#' @exportClass TENxFile
.TENxFile <- setClass(
    Class = "TENxFile",
    contains = "BiocFile",
    slots =
        c(version = "character", group = "character", extension = "character")
)

.check_file_exists <- function(object) {
    if (file.exists(path(object)))
        NULL
    else
        "Path to the file must be valid"
}

.validTENxFile <- function(object) {
    if (length(path(object))) {
        .check_file_exists(object)
    }
}

S4Vectors::setValidity2("TENxFile", .validTENxFile)

#' @exportClass TENxH5
.TENxH5 <- setClass(
    Class = "TENxH5",
    contains = "TENxFile",
    slots = c(extension = "character")
)

.TENxMTX <- setClass(
    Class = "TENxMTX",
    contains = "TENxFile",
    slots = c(extension = "character")
)

.TENxCompressed <- setClass(
    Class = "TENxCompressed",
    contains = "TENxFile",
    slots = c(extension = "character")
)

# TENxFile constructor ----------------------------------------------------

.get_ext <- function(fname) {
    paste0(
        utils::tail(strsplit(basename(fname), "\\.")[[1]], -1),
        collapse = "."
    )
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
#' con <- TENxFile(h5c)
#'
#' @export
TENxFile <- function(
        resource, group = c("matrix", "outs"), version = c("3", "2"), ...
) {
    if (missing(group))
        group <- match.arg(group)
    version <- match.arg(version)
    ext <- .get_ext(resource)
    if (identical(ext, "h5")) {
        l1 <- rhdf5::h5ls(resource, recursive = FALSE)
        gname <- l1[l1$otype == "H5I_GROUP", "name"]
        if (!group %in% gname)
            stop("'group' not found")
    }
    TENxFUN <- switch(
        ext,
        h5 = .TENxH5, mtx = .TENxMTX, tar.gz = .TENxCompressed, .TENxFile
    )
    TENxFUN(
        resource = resource, version = version,
        group = group, extension = ext, ...
    )
}
