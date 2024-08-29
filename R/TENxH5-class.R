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
#'   to version "3". See the internal `data.frame`: `TENxIO:::h5.version.map` for
#'   a map of fields and their corresponding file locations within the H5 file.
#'   This map is used to create the `rowData` structure from the file.
#'
#' @section import:
#'   The `import` method uses `DelayedArray::TENxMatrix` to represent matrix
#'   data. Generally, version 3 datasets contain associated genomic coordinates.
#'   The associated feature data, as displayed by the `rowData` method, is
#'   queried for the "Type" column which will indicate that a `splitAltExps`
#'   operation is appropriate. If a `ref` input is provided to the constructor
#'   function `TENxH5`, it will be used as the main experiment; otherwise, the
#'   most frequent category in the "Type" column will be used. For example,
#'   the Multiome ATAC + Gene Expression feature data contains both 'Gene
#'   Expression' and 'Peaks' labels in the "Type" column.
#'
#' @return A `TENxH5` class object
#'
#' @seealso [TENxH5]
#'
#' @include TENxFile-class.R
#'
#' @exportClass TENxH5
.TENxH5 <- setClass(
    Class = "TENxH5",
    contains = "TENxFile",
    slots = c(version = "character", group = "character", ranges = "character")
)

.get_h5_group <- function(fpath, remote) {
    .checkPkgsAvail("rhdf5")
    l1 <- rhdf5::h5ls(fpath, recursive = FALSE, s3 = remote)
    l1[l1$otype == "H5I_GROUP", "name"]
}

.KNOWN_H5_GROUPS <- c("matrix", "outs")
.KNOWN_VERSIONS <- c("3", "2")

.check_h5_group <- function(group) {
    g_msg <- paste(.KNOWN_H5_GROUPS, collapse = ", ")
    if (!group %in% .KNOWN_H5_GROUPS)
        warning("'group' not in known 10X groups: ", g_msg, call. = FALSE)
}

.getDim <- function(file, group, remote) {
    rhdf5::h5read(file, paste0(group, "/", "shape"), s3 = remote)
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
#'   An additional `ref` argument can be provided when the file contains
#'   multiple `feature_type` in the file or "Type" in the `rowData`. By
#'   default, the most frequent type is represented.
#'
#'   For data that do not contain genomic coordinate information, one can
#'   set the `ranges` argument to `NA_character_`.
#'
#'   The data version "3" mainly includes a "matrix" group and "interval"
#'   information within the file. Version "2" data does not include
#'   ranged-based information and has a different directory structure compared
#'   to version "3". See the internal `data.frame`: `TENxIO:::h5.version.map`
#'   for a map of fields and their corresponding file locations within the H5
#'   file. This map is used to create the `rowData` structure from the file.
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
#'   "/features/interval". Set to `NA_character_` if range information is not
#'   present.
#'
#' @param rowidx,colidx numeric() A vector of indices corresponding to either
#'   rows or columns that will dictate the data imported from the file. The
#'   indices will be passed on to the `[` method of the `TENxMatrix`
#'   representation.
#'
#' @return Usually, a `SingleCellExperiment` instance
#'
#' @seealso `import` section in [TENxH5-class]
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
#' h5f <- system.file(
#'     "extdata", "10k_pbmc_ATACv2_f_bc_ex.h5",
#'     package = "TENxIO", mustWork = TRUE
#' )
#'
#' ## Optional ref input, most frequent Type used by default
#' th5 <- TENxH5(h5f, ranges = "/features/id", ref = "Peaks")
#' th5
#' TENxH5(h5f, ranges = "/features/id")
#' import(th5)
#'
#' @export
TENxH5 <-
    function(resource, version, group, ranges, rowidx, colidx, ...)
{
    remote <- R.utils::isUrl(resource)
    group <- .get_h5_group(resource, remote)
    .check_h5_group(group)
    if (missing(version))
        version <- .get_tenx_version(group)
    dims <- .getDim(resource, group, remote)
    ext <- list(...)[["extension"]]
    if (is.null(ext))
        ext <- .get_ext(resource)
    if (!identical(tolower(ext), "h5"))
        warning("File extension is not 'h5'; import may fail", call. = FALSE)
    if (missing(ranges))
        ranges <- .selectByVersion(h5.version.map, version, "Ranges")
    else if ((!is.character(ranges) && is.na(ranges)) || !nzchar(ranges))
        ranges <- NA_character_
    if (missing(rowidx))
        rowidx <- seq_len(dims[[1L]])
    if (missing(colidx))
        colidx <- seq_len(dims[[2L]])
    .TENxH5(
        resource = resource, group = group, version = version, ranges = ranges,
        rowidx = rowidx,
        colidx = colidx,
        remote = remote,
        extension = ext
    )
}

h5.version.map <- data.frame(
    Version = c("3", "2"),
    ID = c("/features/id", "/genes"),
    Symbol = c("/features/name", "/gene_names"),
    Type = c("/features/feature_type", NA_character_),
    Ranges = c("/features/interval", NA_character_)
)

.selectByVersion <-
    function(df, version, select = !names(df) %in% c("Version", "Ranges"))
{
    if (missing(version) || is.na(version))
        NA_character_
    else
        df[df[["Version"]] == version, select]
}

.read_rowData <- function(x, nrows) {
    gm <- .selectByVersion(h5.version.map, x@version)
    remote <- x@remote
    gm[] <- Filter(Negate(is.na), gm)
    res <- lapply(gm, function(colval) {
        readname <- paste0(x@group, colval)
        as.character(
            rhdf5::h5read(path(x), index = list(nrows), readname, s3 = remote)
        )
    })
    DF <- as(res, "DataFrame")
    if ("Type" %in% names(DF))
        DF[["Type"]] <- as.factor(DF[["Type"]])
    rownames(DF) <- nrows
    DF
}

#' @describeIn TENxH5 Generate the rowData ad hoc from a TENxH5 file
#'
#' @param x A `TENxH5` object
#'
#' @inheritParams SummarizedExperiment::rowData
#'
#' @export
setMethod("rowData", "TENxH5", function(x, use.names = TRUE, ...) {
    nrows <- list(...)[["rows"]]
    ## Implement a smaller index for display purposes only
    mxrow <- max(x@rowidx)
    if (is.null(nrows) && mxrow > 12)
        nrows <- c(seq(6), mxrow - 5:0)
    .read_rowData(x, nrows)
})

#' @describeIn TENxH5 Get the dimensions of the data as stored in the file
#' @export
setMethod("dim", "TENxH5", function(x) {
    c(length(x@rowidx), length(x@colidx))
})

#' @describeIn TENxH5 Get the dimension names from the file
#' @export
setMethod("dimnames", "TENxH5", function(x) {
    id <- .selectByVersion(h5.version.map, x@version, "ID")
    list(
        rhdf5::h5read(path(x), paste0(x@group, "/", id), s3 = x@remote),
        rhdf5::h5read(path(x), paste0(x@group, "/", "barcodes"), s3 = x@remote)
    )
})

#' @describeIn TENxH5 Read genome string from file
#' @importFrom GenomeInfoDb genome genome<-
#' @export
setMethod("genome", "TENxH5", function(x) {
    group <- x@group
    version <- x@version
    remote <- x@remote
    if (is.na(x@ranges))
        stop("'rowRanges' data not available, e.g., in '/features/interval'")
    gens <- rhdf5::h5read(
        path(x), paste0(group, "/", "features/genome"), s3 = remote
    )
    ugens <- unique(gens)
    intervals <- rhdf5::h5read(
        path(x), paste0(group, "/", x@ranges), s3 = remote
    )
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
    remote <- x@remote
    if (is.na(x@ranges))
        stop("'rowRanges' data not available, e.g., in '/features/interval'")
    rows <- list(...)[["rows"]]
    ## Implement a smaller index for display purposes only
    mxrow <- max(x@rowidx)
    if (is.null(rows) && mxrow > 12)
        rows <- c(seq(6), mxrow - 5:0)
    interval <- rhdf5::h5read(
        path(x), paste0(group, x@ranges), list(rows), s3 = remote
    )
    ## Hack to allow NA ranges for later removal (keeping data parallel)
    interval[interval == "NA"] <- "NA_character_:0"
    gr <- as(as.character(interval), "GRanges")
    names(gr) <- rows
    mcols(gr) <- rowData(x, rows = rows)
    genome(gr) <- genome(x)
    gr
})

#' @describeIn TENxH5 Import TENxH5 data as a SingleCellExperiment; see section
#'   below
#'
#' @importFrom MatrixGenerics rowRanges
#' @importFrom BiocBaseUtils isScalarCharacter
#' @importFrom S4Vectors isEmpty
#'
#' @inheritParams BiocIO::import
#'
#' @export
setMethod("import", "TENxH5", function(con, format, text, ...) {
    .checkPkgsAvail("HDF5Array")
    matrixdata <-
        HDF5Array::TENxMatrix(path(con), con@group)[con@rowidx, con@colidx]
    dots <- list(...)
    if (!con@version %in% c("2", "3"))
        stop("Version not supported.")
    sce <- SingleCellExperiment(
        assays = list(counts = matrixdata)
    )
    if (identical(con@version, "3")) {
        if (!is.na(con@ranges)) {
            rr <- rowRanges(con, rows = con@rowidx)
            names(rr) <- mcols(rr)[["ID"]]
            rowRanges(sce) <- rr
            ## remove stand-in NA values
            sce <- sce[seqnames(rr) != "NA_character_", ]
        }
        if (!isEmpty(rowData(con))) {
            rowData(sce) <- .read_rowData(con, con@rowidx)
        }
    }
    types <- rowData(sce)[["Type"]]
    if (!is.null(types)) {
        ref <- dots[["ref"]]
        if (is.null(ref) || is.na(ref))
            ref <- names(which.max(table(types)))
        if (isScalarCharacter(ref))
            sce <- splitAltExps(sce, types, ref = ref)
    }
    sce
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
