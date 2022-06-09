#' @exportClass TENxH5
.TENxH5 <- setClass(
    Class = "TENxH5",
    contains = "TENxFile",
    slots = c(version = "character", group = "character")
)

.get_h5_group <- function(fpath) {
    l1 <- rhdf5::h5ls(fpath, recursive = FALSE)
    l1[l1$otype == "H5I_GROUP", "name"]
}

.KNOWN_H5_GROUPS <- c("matrix", "outs")

.check_h5_group <- function(fpath) {
    .checkPkgsAvail("rhdf5")
    gname <- .get_h5_group(fpath)
    if (!gname %in% .KNOWN_H5_GROUPS)
        stop("'group' not recognized")
    gname
}

.getDim <- function(file) {
    rhdf5::h5read(file, "matrix/shape")
}

#' @rdname TENxH5
#' @title Import H5 files from 10X
#' @aliases TENxH5-class
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
    group <- .check_h5_group(resource)
    dims <- .getDim(resource)
    .TENxH5(
        resource = resource, group = group, version = version,
        rowidx = seq_len(dims[[1L]]),
        colidx = seq_len(dims[[2L]]),
        ...
    )
}

gene.meta <- data.frame(
    version = c("3", "2"),
    ID = c("/features/id", "/genes"),
    Symbol = c("/features/name", "/gene_names"),
    Type = c("/features/feature_type", NA_character_)
)

#' @describeIn TENxH5 Generate the rowData ad hoc from a TENxH5 file
#' @export
setMethod("rowData", "TENxH5", function(x, use.names = TRUE, ...) {
    gm <- gene.meta[
        x@version == gene.meta[["version"]],
        !names(gene.meta) %in% c("group", "version")
    ]
    gm[] <- Filter(Negate(is.na), gm)
    res <- lapply(gm, function(colval) {
        readname <- paste0(x@group, colval)
        rhdf5::h5read(path(x), readname)
    })
    as(res, "DataFrame")
})

#' @describeIn TENxH5 Get the dimensions of the data as stored in the file
#' @export
setMethod("dim", "TENxH5", function(x) {
    c(length(x@rowidx), length(x@colidx))
})

#' @describeIn TENxH5 Get the dimension names from the file
#' @export
setMethod("dimnames", "TENxH5", function(x) {
    gm <- gene.meta[
        x@version == gene.meta[["version"]],
        !names(gene.meta) %in% c("group", "version")
    ]
    IRanges::CharacterList(
        rhdf5::h5read(path(x), file.path(x@group, gm[["ID"]])),
        rhdf5::h5read(path(x), "matrix/barcodes")
    )
})

#' @describeIn TENxH5 Read genome string from file
#' @importFrom GenomeInfoDb genome genome<-
#' @export
setMethod("genome", "TENxH5", function(x) {
    intervals <- rhdf5::h5read(path(x), "matrix/features/interval")
    splitints <- strsplit(intervals, ":", fixed = TRUE)
    seqnames <- vapply(splitints, `[[`, character(1L), 1L)
    if (any(seqnames == "NA"))
        warning("'seqlevels' contain NA values")
    gens <- rhdf5::h5read(path(x), "matrix/features/genome")
    vapply(split(gens, seqnames), unique, character(1L))
})

#' @describeIn TENxH5 Read interval data and represent as GRanges
#' @importFrom S4Vectors mcols<-
#' @export
setMethod("rowRanges", "TENxH5", function(x, ...) {
    interval <- rhdf5::h5read(path(x), "matrix/features/interval")
    interval[interval == "NA"] <- "NA_character_:0"
    gr <- as(as.character(interval), "GRanges")
    mcols(gr) <- rowData(x)
    genome(gr) <- genome(x)
    gr
})

#' @describeIn TENxH5 Import TENxH5 data as a SingleCellExperiment
#' @importFrom MatrixGenerics rowRanges
#' @import SingleCellExperiment
#' @export
setMethod("import", "TENxH5", function(con, format, text, ...) {
    .checkPkgsAvail("HDF5Array")
    matrixdata <- HDF5Array::TENxMatrix(path(con), con@group)
    if (identical(con@version, "3")) {
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = matrixdata), rowRanges = rowRanges(con)
        )
        rownames(sce) <- mcols(sce)[["ID"]]
        ## remove stand-in NA values
        sce <- sce[seqnames(rowRanges(sce)) != "NA_character_", ]
        splitAltExps(sce, rowData(sce)[["Type"]], ref = "Gene Expression")
    } else {
        stop("Version 2 not supported yet.")
    }
})
