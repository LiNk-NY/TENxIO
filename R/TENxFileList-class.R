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
#'
#' @slot listData list() The data in list format
#'
#' @slot extension character() A vector of file extensions for each file
#'
#' @slot compressed logical(1) Whether the file is compressed as `.tar.gz`
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
    features = c("features.tsv.gz", "genes.tsv.gz")
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
#' @param ... A single file path, named arguments corresponding to file paths,
#'   or a list of named file paths
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
#' import(TENxFileList(fl))
#'
#' @export
TENxFileList <- function(..., version, compressed = FALSE) {
    dots <- S4Vectors::SimpleList(...)
    exts <- dots[["extension"]]
    if (length(names(dots)))
        dots <- dots[names(dots) != "extension"]
    undots <- dots[[1L]]
    if (identical(length(dots), 1L)) {
        if (is.character(undots) && is.null(exts))
            exts <- .get_ext(undots)
        if (missing(version) && is.character(undots))
            version <- .version_from_filelist(undots)
        if (is(undots, "TENxFile")) {
            exts <- undots@extension
            version <- undots@version
        }
        if (is.list(undots)) {
            exts <- vapply(undots, .get_ext, character(1L))
            version <- .version_from_fnames(unlist(undots, use.names = FALSE))
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

#' @importFrom utils untar tail
.TENxUntar <- function(con) {
    dir.create(tempdir <- tempfile())
    untar(path(con), exdir = tempdir)
    tempdir
}

#' TSVFile: A class to represent 10x tab separated values files
#'
#' This class is general purpose for reading in tabular data from the
#' 10x Genomics website with the `.tsv` file extension. The class also supports
#' compressed files, i.e., those with the `.tsv.gz` extension.
#'
#' @return A `TSVFile` class object; a `tibble` for the import method
#'
#' @keywords internal
.TSVFile <- setClass(Class = "TSVFile", contains = "TENxFile")

#' @describeIn TSVFile General import function for `tsv` files from 10x;
#'   using `readr::read_tsv` and returning a `tibble` representation
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom readr read_tsv
#'
#' @keywords internal
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
    fnames <- basename(fnames)
    if ("features.tsv.gz" %in% fnames)
        version <- "3"
    else if ("genes.tsv.gz" %in% fnames)
        version <- "2"
    else
        NA_character_
}

.version_from_filelist <- function(tarball) {
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
    if (con@version %in% c("2", "3")) {
        mat <- datalist[["matrix.mtx.gz"]]
        colnames(mat) <- unlist(
            datalist[["barcodes.tsv.gz"]], use.names = FALSE
        )
        feats <- datalist[[features]]
        feats[is.na(feats[["Chr"]]), "Chr"] <- "NA_character_:0"
        rr <- GenomicRanges::makeGRangesFromDataFrame(
            feats, keep.extra.columns = TRUE
        )
        sce <- SingleCellExperiment(
            assays = SimpleList(counts = mat), rowRanges = rr
        )
        sce <- sce[seqnames(rr) != "NA_character_", ]
        splitAltExps(
            sce,
            feats[["Type"]],
            ref = "Gene Expression"
        )
    } else {
        ## return a list for now
        datalist
    }
})
