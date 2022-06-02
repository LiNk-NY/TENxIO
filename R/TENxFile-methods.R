#' @include TENxFileList-class.R
NULL

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

#' Importing peak-annotation files
#'
#' @keywords internal
#'
#' @examples
#'
#' fi <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_atac_peak_annotation.tsv"
#' pa <- .import_peak_anno(fi)
#'
.import_peak_anno <- function(file) {
    stopifnot(endsWith(file, "peak_annotation.tsv"))
    panno <- readr::read_tsv(file, col_types = c("c", "n", "n", "c", "n", "c"))
    makeGRangesFromDataFrame(panno, keep.extra.columns = TRUE)
}

#' Importing fragment files
#'
#' @keywords internal
#'
#' @examples
#'
#' con <- TENxFile("~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5")
#' sce <- import(con)
#' fr <- "~/data/10x/pbmc_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
#' fra <- .import_fragments(fr)
#'
.import_fragments <- function(file, withIndex = TRUE) {
    ## read tabix format with Rsamtools
    stopifnot(endsWith(file, "atac_fragments.tsv.gz"))
    tabindex <- paste0(file, ".tbi")
    if (withIndex)
        stopifnot(file.exists(tabindex))

    tb <- Rsamtools::TabixFile(file, tabindex, yieldSize = NA_integer_)
    # reslist <- scanTabix(tb, param = rowRanges(sce))
    # lapply(Filter(length, reslist), function(x) read.table(textConnection(x)))
    ex <- read.table(textConnection(scanTabix(tb)[[1]]))

    ## names from table in the link below
    ## https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
    names(ex) <- c("chrom", "chromStart", "chromEnd", "barcode", "readSupport")
    ggr <- makeGRangesFromDataFrame(
        ex, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
    )
    ggrl <- as(splitAsList(ggr, mcols(ggr)$barcode), "GRangesList")
    RaggedExperiment(ggrl)
    # compactAssay(erg, i = "readSupport", sparse = TRUE)
}
