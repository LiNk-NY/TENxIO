#' TENxFileList: A list-like representation for TENxFiles
#'
#' @description This class was designed to mainly handle tarballs from
#' 10X Genomics. The typical file extension for these tarballs is `.tar.gz`.
#'
#' @details
#' These tarballs usually contain three files:
#'   1. `matrix.mtx.gz` - the counts matrix
#'   2. `features.tsv.gz` - row metadata usually represented as `rowData`
#'   3. `barcodes.tsv.gz` - column names corresponding to cell barcode
#'   identifiers
#' Note that version '2' includes `genes.tsv.gz` instead of `features.tsv.gz` in
#' version '3'.
#'
#' @slot listData list() The data in list format
#'
#' @slot extension character() A vector of file extensions for each file
#'
#' @slot compressed logical(1) Whether the file is compressed as `.tar.gz`
#'
#' @slot version character(1) The version number of the tarball usually either
#'   '2' or '3'
#'
#' @return A `TENxFileList` class object
#'
#' @exportClass TENxFileList
.TENxFileList <- setClass(
    Class = "TENxFileList",
    contains = "SimpleList",
    slots = c(
        listData = "list",
        extension = "character",
        compressed = "logical",
        version = "character"
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

file.list.map <- data.frame(
    Version = c("3", "2"),
    features = c("features.tsv", "genes.tsv")
)

S4Vectors::setValidity2("TENxFileList", .validTENxFileList)

# Constructor -------------------------------------------------------------

#' TENxFileList: Represent groups of files from 10X Genomic
#'
#' @description This constructor function is meant to handle `.tar.gz` tarball
#'   files from 10X Genomics.
#'
#' @details
#' These tarballs usually contain three files:
#'   1. `matrix.mtx.gz` - the counts matrix
#'   2. `features.tsv.gz` - row metadata usually represented as `rowData`
#'   3. `barcodes.tsv.gz` - column names corresponding to cell barcode
#'   identifiers
#' If all above files are in the tarball, the import method will provide a
#' `SingleCellExperiment`. Otherwise, a simple list of imported data is given.
#' Note that version "3" uses 'features.tsv.gz' and version "2" uses
#' 'genes.tsv.gz'. If known, indicate the `version` argument in the
#' `TENxFileList` constructor function.
#'
#' @param ... Typically, a file path to a tarball archive. Can be named
#'   arguments corresponding to file paths, or a named list of file paths.
#'
#' @param version character(1) The version in the tarball. See details.
#'
#' @param compressed logical(1) Whether or not the file provided is compressed,
#'   usually as `tar.gz` (default FALSE)
#'
#' @return Either a `SingleCellExperiment` or a list of imported data
#'
#' @examples
#'
#' fl <- system.file(
#'     "extdata", "pbmc_granulocyte_sorted_3k_ff_bc_ex_matrix.tar.gz",
#'     package = "TENxIO", mustWork = TRUE
#' )
#'
#' ## Method 1 (tarball)
#' TENxFileList(fl)
#'
#' ## import() method
#' import(TENxFileList(fl))
#'
#'
#' dir.create(tdir <- tempfile())
#' untar(fl, exdir = tdir)
#' files <- list.files(tdir, recursive = TRUE, full.names = TRUE)
#' names(files) <- basename(files)
#' filelist <- lapply(files, TENxFile)
#'
#' ## Method 2 (list of files)
#' TENxFileList(filelist, compressed = FALSE)
#'
#' TENxFileList(as(filelist, "SimpleList"), compressed = FALSE)
#' 
#'
#' ## Method 3 (named arguments)
#' TENxFileList(
#'     barcodes.tsv.gz = TENxFile(files[1]),
#'     features.tsv.gz = TENxFile(files[2]),
#'     matrix.mtx.gz = TENxFile(files[3])
#' )
#'
#' @export
TENxFileList <- function(..., version, compressed = FALSE) {
    dots <- S4Vectors::SimpleList(...)
    exts <- dots[["extension"]]
    if (length(names(dots)))
        dots <- dots[names(dots) != "extension"]
    undots <- dots[[1L]]
    if (identical(length(dots), 1L)) {
        if (is.list(undots) || is(undots, "SimpleList")) {
            exts <- vapply(undots, .get_ext, character(1L))
            version <- .version_from_fnames(names(undots))
            dots <- undots
        } 
        isdir <- try(file.info(undots)[["isdir"]], silent = TRUE) 
        if (inherits(isdir, "try-error") && is(undots, "SimpleList")) {
            exts <- vapply(undots, .get_ext, character(1L))
            version <- .version_from_fnames(names(undots))
        } else if (isdir) {
            undots <- .get_files_from_folder(undots)
            dots <- lapply(undots, TENxFile)
        }
        if (is.character(undots) && is.null(exts))
            exts <- .get_ext(undots)
        if (missing(version) && is.character(undots) &&
                all(endsWith(undots, "tar.gz")))
            version <- .version_from_tarball(undots)
        else if (missing(version) && is.character(undots))
            version <- .version_from_fnames(undots)
        if (is(undots, "TENxFile")) {
            exts <- undots@extension
            version <- undots@version
        }
    } else {
        exts <- vapply(dots, .get_ext, character(1L))
        version <- .version_from_fnames(dots)
    }
    if (identical(exts, "tar.gz"))
        compressed <- TRUE
    .TENxFileList(
        dots, extension = exts, compressed = compressed, version = version
    )
}

.get_files_from_folder <- function(folder) {
    files <- list.files(folder, full.names = TRUE, recursive = TRUE)
    names(files) <- basename(files)
    files
}

#' @importFrom utils untar tail
.TENxUntar <- function(con) {
    dir.create(tempdir <- tempfile())
    untar(path(con), exdir = tempdir)
    tempdir
}

.get_path <- function(object) {
    if (is(object, "TENxFile"))
        path(object)
    else
        object
}

#' @describeIn TENxFileList-class Obtain file paths for all files in the object
#'   as a vector
#'
#' @inheritParams BiocGenerics::path
#'
#' @export
setMethod("path", "TENxFileList", function(object, ...) {
    vapply(object, .get_path, character(1L))
})

#' @describeIn TENxFileList-class An intermediate method for decompressing
#'   (via untar) the contents of a `.tar.gz` file list
#'
#' @importFrom BiocIO decompress
#' @importFrom BiocBaseUtils setSlots
#'
#' @inheritParams BiocIO::decompress
#'
#' @param manager A `ConnectionManager` internal instance; currently not used.
#'
#' @export
setMethod("decompress", "TENxFileList", function(manager, con, ...) {
    res_ext <- con@extension
    if (is.na(res_ext))
        res_ext <- .get_ext(path(con))
    if (con@compressed) {
        if (identical(res_ext, "tar.gz")) {
            tenfolder <- .TENxUntar(con)
            gfolder <- list.files(tenfolder, full.names = TRUE)
            if (file.info(gfolder)$isdir)
                gfiles <- list.files(
                    gfolder, recursive = TRUE, full.names = TRUE
                )
            else
                gfiles <- gfolder
            newlistdata <- lapply(.setNames(gfiles, basename(gfiles)), TENxFile)
            con <- setSlots(
                object = con, listData = newlistdata, compressed = FALSE
            )
        } else {
            stop("Extension type: ", res_ext, " not supported")
        }
    }
    con
})

.version_from_fnames <- function(fnames) {
    if (is(fnames, "SimpleList"))
        fnames <- names(fnames)
    else
        fnames <- basename(fnames)
    if (any(grepl("features.tsv[\\.gz]*", fnames)))
        version <- "3"
    else if (any(grepl("genes.tsv[\\.gz]*", fnames)))
        version <- "2"
    else
        NA_character_
}

.version_from_tarball <- function(tarball) {
    flist <- untar(tarball, list = TRUE)
    .version_from_fnames(flist)
}

#' @describeIn TENxFileList-class Recursively import files within a
#'   `TENxFileList`
#'
#' @inheritParams BiocIO::import
#'
#' @export
setMethod("import", "TENxFileList", function(con, format, text, ...) {
    fdata <- decompress(con = con)
    datalist <- lapply(fdata, import)
    features <-
        .selectByVersion(file.list.map, version = con@version, "features")
    features <- grep(features, names(fdata), fixed = TRUE, value = TRUE)
    matrix <- grep("matrix.mtx", names(fdata), fixed = TRUE, value = TRUE)
    barcodes <- grep("barcodes.tsv", names(fdata), fixed = TRUE, value = TRUE)
    if (con@version %in% c("2", "3")) {
        mat <- datalist[[matrix]]
        colnames(mat) <- unlist(
            datalist[[barcodes]], use.names = FALSE
        )
        feats <- datalist[[features]]
        sce <- as(mat, "SingleCellExperiment")
        if ("Chr" %in% names(feats)) {
            feats[is.na(feats[["Chr"]]), "Chr"] <- "NA_character_:0"
            rr <- GenomicRanges::makeGRangesFromDataFrame(
                feats, keep.extra.columns = TRUE
            )
            rowRanges(sce) <- rr
            sce <- sce[seqnames(rr) != "NA_character_", ]
        } else {
            rowData(sce) <- feats
        }
        if (!is.null(feats[["Type"]]) && length(feats[["Type"]]))
            splitAltExps(
                sce,
                feats[["Type"]],
                ref = "Gene Expression"
            )
        else
            sce
    } else {
        ## return a list for now
        datalist
    }
})
