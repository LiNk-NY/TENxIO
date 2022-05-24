#' @importFrom MatrixGenerics rowRanges
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
    gen <- unique(rhdf5::h5read(h5f, gloc))
    if (length(gen) != 1L)
        stop("The genome build in ", gloc, " is not consistent")
    gen
}

#' @importFrom GenomeInfoDb genome genome<-
.getRowRanges <- function(con) {
    interval <- rhdf5::h5read(path(con), "matrix/features/interval")
    interval[interval == "NA"] <- "NA_character_:0"
    gr <- as(as.character(interval), "GRanges")
    mcols(gr) <- .getRowDat(con)
    genbuild <- rep(.getGenome(con), length(genome(gr)))
    genome(gr) <- genbuild
    gr
}

setMethod("import", "TENxH5", function(con, format, text, ...) {
    matrixdata <- HDF5Array::TENxMatrix(path(con), con@group)
    if (!requireNamespace("rhdf5", quietly = TRUE))
        stop("Install 'rhdf5' to import TENxH5 data")
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrixdata), rowRanges = .getRowRanges(con)
    )
})

setMethod("import", "TENxMTX", function(con, format, text, ...) {
    mtxf <- SingleCellMultiModal:::.read_mtx(path(con))
    ## TODO: make use of other files
    SummarizedExperiment(assays = SimpleList(counts = mtxf))
})


.TENxUntar <- function(con) {
    dir.create(tempdir <- tempfile())
    untar(con, exdir = tempdir)
    tempdir
}

.TENxDecompress <- function(con) {
    res_ext <- .get_ext(resource)
    if (identical(res_ext, "tar.gz")) {
        tenfolder <- .TENxUntar(con)
        gfiles <- list.files(tenfolder, full.names = TRUE)
        ## Sort through files and import
        ## import()
    } else {
        stop("Extension type: ", res_ext, " not supported")
    }
}

setMethod("import", "TENxCompressed", function(con, format, text, ...) {
   ## Decompress and callNextMethod to proper file type
})
