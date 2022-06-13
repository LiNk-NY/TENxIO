#' TENxH5: The HDF5 file representation class for 10X Data
#'
#' @description This class is designed to work with 10x Single Cell datasets.
#'   It was developed using the PBMC 3k 10X dataset from the CellRanger v2
#'   pipeline.
#'
#' @slot version character(1) There are currently two recognized versions
#'   associated with 10X data, either version "2" or "3". See details for more
#'   information.
#'
#' @slot group character(1) The HDF5 group embedded within the file structure,
#'   this is usually either the "matrix" or "outs" group but other groups are
#'   supported as well.
#'
#' @details The data version "3" mainly includes a "matrix" group and "interval"
#'   information within the file. Version "2" data does not include
#'   ranged-based information and has a different directory structure compared
#'   to version "3". See the internal `data.frame`: `TENxIO:::gene.meta` for
#'   a map of fields and their corresponding file locations within the H5 file.
#'   This map is used to create the `rowData` structure from the file.
#'
#' @include TENxFile-class.R
#' @exportClass TENxH5
.TENxH5 <- setClass(
    Class = "TENxH5",
    contains = "TENxFile",
    slots = c(version = "character", group = "character")
)

.get_h5_group <- function(fpath) {
    .checkPkgsAvail("rhdf5")
    l1 <- rhdf5::h5ls(fpath, recursive = FALSE)
    l1[l1$otype == "H5I_GROUP", "name"]
}

.KNOWN_H5_GROUPS <- c("matrix", "outs")
.KNOWN_VERSIONS <- c("3", "2")

.check_h5_group <- function(group) {
    g_msg <- paste(.KNOWN_H5_GROUPS, collapse = ", ")
    if (!group %in% .KNOWN_H5_GROUPS)
        warning("'group' not in known 10X groups: ", g_msg)
}

.getDim <- function(file, group) {
    rhdf5::h5read(file, paste0(group, "/", "shape"))
}

.get_tenx_version <- function(group) {
    .KNOWN_VERSIONS[match(group, .KNOWN_H5_GROUPS)]
}

#' TENxH5: Import H5 files from 10X
#'
#' This constructor function was developed using the PBMC 3K dataset from 10X
#' Genomics (verison 3). Other versions are supported and input arguments
#' `version` and `group` can be overridden.
#'
#' @inheritParams TENxFile
#'
#' @param version character(1) There are currently two recognized versions
#'   associated with 10X data, either version "2" or "3". See details for more
#'   information.
#'
#' @param group character(1) The HDF5 group embedded within the file structure,
#'   this is usually either the "matrix" or "outs" group but other groups are
#'   supported as well (e.g., "mm10").
#'
#' @details The data version "3" mainly includes a "matrix" group and "interval"
#'   information within the file. Version "2" data does not include
#'   ranged-based information and has a different directory structure compared
#'   to version "3". See the internal `data.frame`: `TENxIO:::gene.meta` for
#'   a map of fields and their corresponding file locations within the H5 file.
#'   This map is used to create the `rowData` structure from the file.
#'
#' @examples
#'
#' h5f <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
#' con <- TENxH5(h5f)
#' import(con)
#'
#' @export
TENxH5 <-
    function(resource, version, group, ...)
{
    group <- .get_h5_group(resource)
    .check_h5_group(group)
    if (missing(version))
        version <- .get_tenx_version(group)
    dims <- .getDim(resource, group)
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
        as.character(rhdf5::h5read(path(x), readname))
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
    list(
        rhdf5::h5read(path(x), paste0(x@group, "/", gm[["ID"]])),
        rhdf5::h5read(path(x), paste0(x@group, "/", "barcodes"))
    )
})

#' @describeIn TENxH5 Read genome string from file
#' @importFrom GenomeInfoDb genome genome<-
#' @export
setMethod("genome", "TENxH5", function(x) {
    group <- x@group
    gens <- rhdf5::h5read(path(x), paste0(group, "/", "features/genome"))
    ugens <- unique(gens)
    intervals <- rhdf5::h5read(path(x), paste0(group, "/", "features/interval"))
    splitints <- strsplit(intervals, ":", fixed = TRUE)
    seqnames <- vapply(splitints, `[[`, character(1L), 1L)
    if (any(seqnames == "NA"))
        warning("'seqlevels' contain NA values")
    if (identical(length(ugens), 1L)) {
        useq <- unique(seqnames)
        .setNames(rep(ugens, length(useq)), useq)
    } else
        vapply(split(gens, seqnames), unique, character(1L))
})

#' @describeIn TENxH5 Read interval data and represent as GRanges
#' @importFrom S4Vectors mcols<-
#' @export
setMethod("rowRanges", "TENxH5", function(x, ...) {
    group <- x@group
    interval <- rhdf5::h5read(path(x), paste0(group, "/features/interval"))
    ## Hack to allow NA ranges for later removal (keeping data parallel)
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
