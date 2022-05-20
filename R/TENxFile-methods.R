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

.getRowRanges <- function(con) {
    interval <- rhdf5::h5read(path(con), "matrix/features/interval")
    interval[interval == "NA"] <- "NA_character_:0"
    gr <- as(as.character(interval), "GRanges")
    mcols(gr) <- .getRowDat(con)
    gr
}

setMethod("import", "TENxFile", function(con, format, text, ...) {
    matrixdata <- HDF5Array::TENxMatrix(path(con), con@group)
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrixdata), rowRanges = .getRowRanges(con)
    )
})
