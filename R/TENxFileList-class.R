#' TENxFileList: A list-like representation for TENxFiles
#'
#' @description This class was designed to mainly handle tarballs from
#' 10X Genomics. The typical file extension for these tarballs is `.tar.gz`.
#'
#' @details
#' These tarballs usually contain three files:
#'   1. `matrix.mtx.gz` - the counts matrix
#'   2. `features.tsv.gz` - row metadata usually represented as `rowData`
#'   3. `barcodes.tsv.gz` - column names corresponding to cell barcode identifiers
#'
#' @slot listData list() The data in list format
#' @slot extension character() A vector of file extensions for each file
#' @slot compressed logical(1) Whether the file is compressed as `.tar.gz`
#'
#' @exportClass TENxFileList
.TENxFileList <- setClass(
    Class = "TENxFileList",
    contains = "SimpleList",
    slots = c(
        listData = "list", extension = "character", compressed = "logical"
    )
)

.validTENxFileList <- function(object) {
    if (!object@compressed) {
        validClasses <- vapply(object, function(ifile) {
            is(ifile, "TENxFile")
        }, logical(1L))
        if (all(validClasses))
            TRUE
        else
            "Some files are not of class 'TENxFile'"
    }
}

S4Vectors::setValidity2("TENxFileList", .validTENxFileList)

#' @examples
#'
#' fl <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.tar.gz"
#' con <- TENxFileList(fl)
#'
#' @export
TENxFileList <- function(..., compressed = FALSE) {
    dots <- S4Vectors::SimpleList(...)
    undots <- dots[[1L]]
    if (identical(length(dots), 1L)) {
        if (is.character(undots))
            exts <- .get_ext(undots)
        if (is(undots, "TENxFile"))
            exts <- undots@extension
        if (is.list(undots))
            exts <- vapply(undots, .get_ext, character(1L))
    } else {
        exts <- vapply(dots, .get_ext, character(1L))
    }
    if (identical(exts, "tar.gz"))
        compressed <- TRUE
    .TENxFileList(dots, extension = exts, compressed = compressed)
}

#' @importFrom utils untar tail
.TENxUntar <- function(con) {
    dir.create(tempdir <- tempfile())
    untar(path(con), exdir = tempdir)
    tempdir
}

.TSVFile <- setClass(Class = "TSVFile", contains = "TENxFile")

setMethod("import", "TSVFile", function(con, format, text, ...) {
    resource <- path(con)
    df <- readr::read_tsv(
        resource, col_names = FALSE, show_col_types = FALSE, progress = FALSE,
        ...
    )
    fname <- basename(resource)
    if (identical(fname, "features.tsv.gz"))
        names(df) <- c("ID", "Symbol", "Type", "Chr", "Start", "End")
    else if (identical(fname, "barcodes.tsv.gz"))
        names(df) <- "barcode"
    df
})

.get_path <- function(object) {
    if (is(object, "TENxFile"))
        path(object)
    else
        object
}

setMethod("path", "TENxFileList", function(object, ...) {
    vapply(object, .get_path, character(1L))
})

.TENxDecompress <- function(con) {
    res_ext <- .get_ext(path(con))
    if (!con@compressed)
        stop("<internal> ", class(con), " not compressed")
    if (identical(res_ext, "tar.gz")) {
        tenfolder <- .TENxUntar(con)
        gfolder <- list.files(tenfolder, full.names = TRUE)
        if (file.info(gfolder)$isdir)
            gfiles <- list.files(gfolder, recursive = TRUE, full.names = TRUE)
        else
            gfiles <- gfolder
        lapply(gfiles, TENxFile)
    } else {
        stop("Extension type: ", res_ext, " not supported")
    }
}

.TARFILENAMES <- c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")

#' @examples
#'
#' fl <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.tar.gz"
#' con <- TENxFileList(fl)
#' import(con)
#'
#' @export
setMethod("import", "TENxFileList", function(con, format, text, ...) {
    if (con@compressed)
        fdata <- .TENxDecompress(con)
    else
        fdata <- con@listData
    fldata <- TENxFileList(fdata)
    names(fdata) <- basename(path(fldata))
    datalist <- lapply(fdata, import)
    if (all(.TARFILENAMES %in% names(datalist))) {
        mat <- datalist[["matrix.mtx.gz"]]
        colnames(mat) <- unlist(
            datalist[["barcodes.tsv.gz"]], use.names = FALSE
        )
        warning("Matrix of mixed types; see in rowData(x)$Type")
        SingleCellExperiment::SingleCellExperiment(
            SimpleList(counts = mat),
            rowData = datalist[["features.tsv.gz"]]
        )
    } else {
        ## return a list for now
        datalist
    }
})
