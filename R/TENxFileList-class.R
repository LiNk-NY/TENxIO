#' @exportClass TENxFileList
.TENxFileList <- setClass(
    Class = "TENxFileList",
    contains = "SimpleList",
    slots = c(compressed = "logical")
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

#' @export
TENxFileList <- function(..., compressed = FALSE) {
    dots <- list(...)
    if (identical(length(dots[[1]]), 1L) && is.character(dots[[1]]))
        exts <- .get_ext(dots[[1]])
    else if (is.list(dots[[1]]))
        exts <- vapply(dots[[1]], .get_ext, character(1L))
    ## TODO: handle folder and list of files, in addition to tar.gz
    if (identical(exts, "tar.gz"))
        compressed <- TRUE
    asl <- SimpleList(dots[[1]])
    .TENxFileList(asl, compressed = compressed)
}

#' @importFrom utils untar tail
.TENxUntar <- function(con) {
    dir.create(tempdir <- tempfile())
    untar(path(con), exdir = tempdir)
    tempdir
}

.readInFuns <- function(files) {
    file_exts <- .get_ext(files)
    lapply(.setNames(file_exts, basename(files)), function(ext) {
        switch(
            ext,
            mtx.gz = Matrix::readMM,
            tsv.gz = function(...)
                readr::read_tsv(col_names = FALSE, show_col_types = FALSE, ...)
        )
    })
}

.cleanUpFuns <- function(datalist) {
    if (is.null(names(datalist)))
        stop("'datalist' names must correspond to originating file names")
    lapply(.setNames(nm = names(datalist)), function(fname) {
        switch(
            fname,
            features.tsv.gz = function(df) {
                names(df) <- c("ID", "Symbol", "Type", "Chr", "Start", "End")
                df
            },
            barcodes.tsv.gz = function(df) {
                names(df) <- "barcode"
                df
            },
            matrix.mtx.gz = function(mat) {
                as(mat, "dgCMatrix")
            },
        )
    })
}

.TENxDecompress <- function(con) {
    res_ext <- .get_ext(path(con))
    if (identical(res_ext, "tar.gz")) {
        tenfolder <- .TENxUntar(con)
        gfolder <- list.files(tenfolder, full.names = TRUE)
        if (file.info(gfolder)$isdir)
            gfiles <- list.files(gfolder, recursive = TRUE, full.names = TRUE)
        else
            gfiles <- gfolder
        gdata <- Map(
            f = function(reader, x) {
                reader(x)
            }, reader = .readInFuns(gfiles), x = gfiles
        )
        Map(f = function(cleaner, x) {
                cleaner(x)
            }, cleaner = .cleanUpFuns(gdata), x = gdata
        )
    } else {
        stop("Extension type: ", res_ext, " not supported")
    }
}

#' @export
setMethod("import", "TENxFileList", function(con, format, text, ...) {
    if (con@compressed)
        fdata <- .TENxDecompress(con)
    else
        fdata <- con@listData
    mat <- fdata[["matrix.mtx.gz"]]
    colnames(mat) <- unlist(fdata[["barcodes.tsv.gz"]])
    warning("Matrix of mixed types; see in rowData(x)")
    SingleCellExperiment::SingleCellExperiment(
        SimpleList(counts = mat),
        rowData = fdata[["features.tsv.gz"]]
    )
})
