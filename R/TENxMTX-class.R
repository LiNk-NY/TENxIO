.TENxMTX <- setClass(
    Class = "TENxMTX",
    contains = "TENxFile"
)

#' @import SummarizedExperiment
#' @export
setMethod("import", "TENxMTX", function(con, format, text, ...) {
    mtxf <- SingleCellMultiModal:::.read_mtx(path(con))
    ## TODO: make use of other files
    SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(counts = mtxf)
    )
})
