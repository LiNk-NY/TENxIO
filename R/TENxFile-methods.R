#' @importFrom MatrixGenerics rowRanges
#' @include TENxFileList-class.R
NULL

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
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrixdata), rowRanges = .getRowRanges(con)
    )
})

#' @import SummarizedExperiment
#' @export
setMethod("import", "TENxMTX", function(con, format, text, ...) {
    mtxf <- SingleCellMultiModal:::.read_mtx(path(con))
    ## TODO: make use of other files
    SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(counts = mtxf)
    )
})


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
