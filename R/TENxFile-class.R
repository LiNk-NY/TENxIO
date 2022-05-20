#' @importClassesFrom BiocIO BiocFile
#' @importFrom BiocIO import
#'
#' @exportClass TENxFile
.TENxFile <- setClass(
    Class = "TENxFile",
    contains = "BiocFile",
    slots = c(version = "character", group = "character")
)

# TENxFile constructor ----------------------------------------------------

#' TENxFile constructor function
#'
#' @examples
#'
#' h5f <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
#' con <- TENxFile(h5f)
#' import(con)
#'
#' @export
TENxFile <- function(
        resource, group = c("matrix", "outs"), version = c("3", "2")
) {
    if (missing(group))
        group <- match.arg(group)
    version <- match.arg(version)
    ext <- tools::file_ext(resource)
    l1 <- rhdf5::h5ls(resource, recursive = FALSE)
    gname <- l1[l1$otype == "H5I_GROUP", "name"]
    if (!group %in% gname)
        stop("'group' not found")
    if (identical(ext, "h5"))
        .TENxH5(
            resource = resource, version = version,
            group = group, extension = ext
        )
    else
        .TENxFile(
            resource = resource, version = version, group = group
        )
}

#' @exportClass TENxH5
.TENxH5 <- setClass(
    Class = "TENxH5", contains = "TENxFile", slots = c(extension = "character")
)
