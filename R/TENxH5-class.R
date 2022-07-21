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
#' @slot ranges character(1) The HDF5 internal folder location embedded within
#'   the file that points to the ranged data information, e.g.,
#'   "/features/interval".
#'
#' @details The data version "3" mainly includes a "matrix" group and "interval"
#'   information within the file. Version "2" data does not include
#'   ranged-based information and has a different directory structure compared
#'   to version "3". See the internal `data.frame`: `TENxIO:::gene.meta` for
#'   a map of fields and their corresponding file locations within the H5 file.
#'   This map is used to create the `rowData` structure from the file.
#'
#' @include TENxFile-class.R
#'
#' @exportClass TENxH5
.TENxH5 <- setClass(
    Class = "TENxH5",
    contains = "TENxFile",
    slots = c(version = "character", group = "character", ranges = "character")
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
        warning("'group' not in known 10X groups: ", g_msg, call. = FALSE)
}

.getDim <- function(file, group) {
    rhdf5::h5read(file, paste0(group, "/", "shape"))
}

.get_tenx_version <- function(group) {
    .KNOWN_VERSIONS[match(group, .KNOWN_H5_GROUPS)]
}

# Constructor -------------------------------------------------------------

#' TENxH5: Represent H5 files from 10X
#'
#' This constructor function was developed using the PBMC 3K dataset from 10X
#' Genomics (version 3). Other versions are supported and input arguments
#' `version` and `group` can be overridden.
#'
#' @details The various `TENxH5` methods including `rowData` and `rowRanges`,
#'   provide a snapshot of the data using a length 12 head and tail subset for
#'   efficiency. In contrast, methods such as `dimnames` and `dim` give
#'   a full view of the dimensions of the data. The `show` method provides
#'   relevant information regarding the dimensions of the data including
#'   metadata such as `rowData` and "Type" column, if available. The term
#'   "projection" refers to the data class that will be provided once the
#'   data file is `import`ed.
#'
#'   The data version "3" mainly includes a "matrix" group and "interval"
#'   information within the file. Version "2" data does not include
#'   ranged-based information and has a different directory structure compared
#'   to version "3". See the internal `data.frame`: `TENxIO:::gene.meta` for
#'   a map of fields and their corresponding file locations within the H5 file.
#'   This map is used to create the `rowData` structure from the file.
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
#' @param ranges character(1) The HDF5 internal folder location embedded within
#'   the file that points to the ranged data information, e.g.,
#'   "/features/interval".
#'
#' @examples
#'
#' h5f <- system.file(
#'     "extdata", "pbmc_granulocyte_ff_bc_ex.h5",
#'     package = "TENxIO", mustWork = TRUE
#' )
#'
#' TENxH5(h5f)
#'
#' import(TENxH5(h5f))
#'
#' h52 <- "~/data/10x/pbmc_10k/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5"
#' con <- TENxH5(h52, ranges = "/features/id")
#'
#' @export
TENxH5 <-
    function(resource, version, group, ranges, ...)
{
    group <- .get_h5_group(resource)
    .check_h5_group(group)
    if (missing(version))
        version <- .get_tenx_version(group)
    dims <- .getDim(resource, group)
    ext <- list(...)[["extension"]]
    if (is.null(ext))
        ext <- .get_ext(resource)
    if (!identical(tolower(ext), "h5"))
        warning("File extension is not 'h5'; import may fail", call. = FALSE)
    if (missing(ranges))
        ranges <- .selectByVersion(gene.meta, version, "Ranges")
    .TENxH5(
        resource = resource, group = group, version = version, ranges = ranges,
        rowidx = seq_len(dims[[1L]]),
        colidx = seq_len(dims[[2L]]),
        extension = ext
    )
}

gene.meta <- data.frame(
    Version = c("3", "2"),
    ID = c("/features/id", "/genes"),
    Symbol = c("/features/name", "/gene_names"),
    Type = c("/features/feature_type", NA_character_),
    Ranges = c("/features/interval", NA_character_)
)

.selectByVersion <-
    function(df, version, select = !names(df) %in% c("Version", "Ranges"))
{
    df[df[["Version"]] == version, select]
}

#' @describeIn TENxH5 Generate the rowData ad hoc from a TENxH5 file
#'
#' @inheritParams SummarizedExperiment::rowData
#'
#' @export
setMethod("rowData", "TENxH5", function(x, use.names = TRUE, ...) {
    gm <- .selectByVersion(gene.meta, x@version)
    nrows <- list(...)[["rows"]]
    ## Implement a smaller index for display purposes only
    mxrow <- max(x@rowidx)
    if (is.null(nrows) && mxrow > 12)
        nrows <- c(1:6, mxrow - 5:0)
    gm[] <- Filter(Negate(is.na), gm)
    res <- lapply(gm, function(colval) {
        readname <- paste0(x@group, colval)
        as.character(rhdf5::h5read(path(x), index = list(nrows), readname))
    })
    DF <- as(res, "DataFrame")
    if ("Type" %in% names(DF))
        DF[["Type"]] <- as.factor(DF[["Type"]])
    rownames(DF) <- nrows
    DF
})

#' @describeIn TENxH5 Get the dimensions of the data as stored in the file
#' @export
setMethod("dim", "TENxH5", function(x) {
    c(length(x@rowidx), length(x@colidx))
})

#' @describeIn TENxH5 Get the dimension names from the file
#' @export
setMethod("dimnames", "TENxH5", function(x) {
    id <- .selectByVersion(gene.meta, x@version, "ID")
    list(
        rhdf5::h5read(path(x), paste0(x@group, "/", id)),
        rhdf5::h5read(path(x), paste0(x@group, "/", "barcodes"))
    )
})

#' @describeIn TENxH5 Read genome string from file
#' @importFrom GenomeInfoDb genome genome<-
#' @export
setMethod("genome", "TENxH5", function(x) {
    group <- x@group
    version <- x@version
    if (is.na(x@ranges))
        stop("'rowRanges' data not available, e.g., in '/features/interval'")
    gens <- rhdf5::h5read(path(x), paste0(group, "/", "features/genome"))
    ugens <- unique(gens)
    intervals <- rhdf5::h5read(path(x), paste0(group, "/", x@ranges))
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
    version <- x@version
    if (is.na(x@ranges))
        stop("'rowRanges' data not available, e.g., in '/features/interval'")
    rows <- list(...)[["rows"]]
    ## Implement a smaller index for display purposes only
    mxrow <- max(x@rowidx)
    if (is.null(rows) && mxrow > 12)
        rows <- c(1:6, mxrow - 5:0)
    interval <- rhdf5::h5read(
        path(x), paste0(group, x@ranges), list(rows)
    )
    ## Hack to allow NA ranges for later removal (keeping data parallel)
    interval[interval == "NA"] <- "NA_character_:0"
    gr <- as(as.character(interval), "GRanges")
    names(gr) <- rows
    mcols(gr) <- rowData(x, rows = rows)
    genome(gr) <- genome(x)
    gr
})

#' @describeIn TENxH5 Import TENxH5 data as a SingleCellExperiment
#'
#' @importFrom MatrixGenerics rowRanges
#'
#' @inheritParams BiocIO::import
#'
#' @export
setMethod("import", "TENxH5", function(con, format, text, ...) {
    .checkPkgsAvail("HDF5Array")
    matrixdata <- HDF5Array::TENxMatrix(path(con), con@group)
    if (identical(con@version, "3")) {
        rr <- rowRanges(con, rows = con@rowidx)
        sce <- SingleCellExperiment(
            assays = list(counts = matrixdata),
            rowRanges = rr
        )
        rownames(sce) <- mcols(sce)[["ID"]]
        ## remove stand-in NA values
        sce <- sce[seqnames(rr) != "NA_character_", ]
        splitAltExps(
            sce,
            rowData(sce, rows = con@rowidx)[["Type"]],
            ref = "Gene Expression"
        )
    } else {
        stop("Version 2 not supported yet.")
    }
})

#' @describeIn TENxH5 Display a snapshot of the contents within a TENxH5 file
#'   before import
#'
#' @param object A `TENxH5` class object
#'
#' @importFrom BiocBaseUtils selectSome
#' @export
setMethod("show", "TENxH5", function(object) {
    rn <- cn <- NULL
    rno <- rownames(object)
    cno <- colnames(object)
    if (length(rno))
        rn <- selectSome(rno)
    if (length(cno))
        cn <- selectSome(cno)
    rd <- rowData(object)
    rdnames <- selectSome(names(rd))
    rdnamecount <- paste0("names(", length(rd), "):")

    cat(
        class(object), "object",
        "\nresource:", path(object),
        "\nprojection:", "SingleCellExperiment",
        "\ndim:", dim(object),
        "\nrownames:", rn,
        "\nrowData", rdnamecount, rdnames,
        if (length(rd[["Type"]]))
            "\n  Type:", selectSome(levels(rd[["Type"]])),
        "\ncolnames:", cn,
        "\n", sep = " "
    )
})
