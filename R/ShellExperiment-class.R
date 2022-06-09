#' @include TENxH5-class.R
#' @exportClass ShellExperiment
.ShellExperiment <- setClass(
    Class = "ShellExperiment",
    contains = "TENxH5",
    slots = c(
        assays = "Assays_OR_NULL", NAMES = "character_OR_NULL",
        rowRanges = "GenomicRanges_OR_GRangesList", colData = "DataFrame",
        elementMetadata = "DataFrame"
    )
)

#' @export
setMethod("show", "ShellExperiment", function(object) {
    rn <- cn <- NULL
    rno <- rownames(object)
    cno <- colnames(object)
    if (length(rno))
        rn <- BiocBaseUtils::selectSome(rno)
    if (length(cno))
        cn <- BiocBaseUtils::selectSome(cno)
    rd <- rowData(object)
    rdn <- BiocBaseUtils::selectSome(names(rd))

    cat(
        "class:", class(object),
        "\nprojection:", "SingleCellExperiment",
        "\ndim:", dim(object),
        "\nrownames:", rn,
        paste0("\nrowData names(", length(rd), "):"), rdn,
        "\ncolnames:", cn,
        "\n", sep = " "
    )
})

#' @examples
#'
#' h5f <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
#' con <- TENxFile(h5f)
#' ShellExperiment(con)
#'
#' @export
ShellExperiment <- function(con) {
    .ShellExperiment(
        con,
        assays = NULL,
        rowRanges = rowRanges(con),
        elementMetadata = DataFrame(),
        colData = DataFrame()
    )
}
