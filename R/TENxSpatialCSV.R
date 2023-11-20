.TENxSpatialCSV <- setClass(
    Class = "TENxSpatialCSV",
    contains = "TENxFile",
    slots = c(
        isList = "logical",
        colnames = "character"
    )
)

.TISSUE_POS_COLS <- c(
    "barcode", "in_tissue", "array_row", "array_col",
    "pxl_row_in_fullres", "pxl_col_in_fullres"
)

TENxSpatialCSV <- function(resource, colnames = .TISSUE_POS_COLS) {
    if (!is(resource, "TENxFile"))
        resource <- TENxFile(resource)
    isList <- grepl("_list", path(resource), fixed = TRUE)
    .TENxSpatialCSV(
        resource, isList = isList, colnames = colnames
    )
}

setMethod("import", "TENxSpatialCSV", function(con, format, text, ...) {
    read.csv(
        path(con),
        header = !con@isList,
        row.names = 1L,
        col.names = con@colnames
    )
})
