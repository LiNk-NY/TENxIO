#' The ShellExperiment class for temporary representations of H5 files
#'
#' @description The class represents a deconstructed `SingleCellExperiment`
#'   class for inspecting several aspects of the `TENxH5` file class. It
#'   inherits the methods from the `TENxH5` class. See the `TENxH5-class`
#'   documentation.
#'
#' @slot assays Assays_OR_NULL; Assay class input similar to that in the
#'   `SummarizedExperiment` constructor (default `NULL`)
#'
#' @slot rowRanges GenomicRanges_OR_GRangesList; a ranged based representation
#'   of the genomic data within the file
#'
#' @include TENxH5-class.R
#' @exportClass ShellExperiment
.ShellExperiment <- setClass(
    Class = "ShellExperiment",
    contains = "TENxH5",
    slots = c(
        assays = "Assays_OR_NULL",
        rowRanges = "GenomicRanges_OR_GRangesList"
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

#' ShellExperiment constructor function
#'
#' @description The `ShellExperiment` function takes a version 2 or 3 HDF5
#'   dataset typically obtained via the 10x website and provides an intermediate
#'   representation of the data within the file. Supported methods include
#'   `rowRanges`, `dimnames`, `dim`, `rowData`, and `genome`. Currently,
#'   only data with interval information is supported (version 3).
#'
#' @param con character(1) A path to an H5 file
#'
#' @examples
#'
#' h5f <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
#' con <- TENxFile(h5f)
#' ShellExperiment(con)
#'
#' ## Alternatively
#' mat <- HDF5Array::TENxMatrix(h5f, "matrix")
#' hub <- ExperimentHub::ExperimentHub()
#' fname <- hub[["EH1039"]]
#' x5 <- TENxH5(fname, group = "mm10", version = "2")
#'
#' @export
ShellExperiment <- function(con) {
    .ShellExperiment(
        con,
        assays = NULL,
        rowRanges = rowRanges(con)
    )
}
