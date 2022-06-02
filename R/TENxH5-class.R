#' @importFrom MatrixGenerics rowRanges
NULL

#' @exportClass TENxH5
.TENxH5 <- setClass(
    Class = "TENxH5",
    contains = "TENxFile",
    slots = c(version = "character", group = "character")
)

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

gene.meta <- data.frame(
    version = c("3", "2"),
    ID = c("/features/id", "/genes"),
    Symbol = c("/features/name", "/gene_names"),
    Type = c("/features/feature_type", NA_character_)
)

.getRowDat <- function(con) {
    gm <- gene.meta[
        con@version == gene.meta[["version"]],
        !names(gene.meta) %in% c("group", "version")
    ]
    gm[] <- Filter(Negate(is.na), gm)
    res <- lapply(gm, function(colval) {
        readname <- paste0(con@group, colval)
        rhdf5::h5read(path(con), readname)
    })
    as.data.frame(res)
}

.getGenome <- function(con) {
    gloc <- "matrix/features/genome"
    gen <- unique(rhdf5::h5read(path(con), gloc))
    if (length(gen) != 1L)
        stop("The genome build in ", gloc, " is not consistent")
    gen
}

#' @importFrom GenomeInfoDb genome genome<-
#' @importFrom S4Vectors mcols<-
.getRowRanges <- function(con) {
    interval <- rhdf5::h5read(path(con), "matrix/features/interval")
    interval[interval == "NA"] <- "NA_character_:0"
    gr <- as(as.character(interval), "GRanges")
    mcols(gr) <- .getRowDat(con)
    genbuild <- rep(.getGenome(con), length(genome(gr)))
    genome(gr) <- genbuild
    gr
}

#' @import SingleCellExperiment
#' @export
setMethod("import", "TENxH5", function(con, format, text, ...) {
    if (!requireNamespace("HDF5Array", quietly = TRUE))
        stop("Install 'HDF5Array' to import TENxH5 files")
    matrixdata <- HDF5Array::TENxMatrix(path(con), con@group)
    if (identical(con@version, "3")) {
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = matrixdata), rowRanges = .getRowRanges(con)
        )
        rownames(sce) <- mcols(sce)[["ID"]]
        splitAltExps(sce, rowData(sce)[["Type"]], ref = "Gene Expression")
    } else {
        stop("Version 2 not supported yet.")
    }
})
