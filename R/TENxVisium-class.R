.TENxVisium <- setClass(
    Class = "TENxVisium",
    slots = c(
        resources = "TENxFileList",
        spatialList = "TENxSpatialList",
        coordNames = "character"
    )
)

TENxVisium <- function(
    resources,
    spatialResource,
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres")
) {
    if (!is(resources, "TENxFileList"))
        resources <- TENxFileList(resources)
    if (!is(spatialResource, "TENxSpatialList"))
        spatialResource <- TENxSpatialList(spatialResource)

    .TENxVisium(
        resources = resources,
        spatialList = spatialResource,
        coordNames = spatialCoordsNames
    )
}


.validTENxVisium <- function(object) {
    isFL <- is(object@resources, "TENxFileList")
    if (!isFL)
        "'TENxFileList' component is not of TENxFileList class"
    else
        isFL
    isSL <- is(object@spatialList, "TENxSpatialList")
    if (!isSL)
        "'TENxSpatialList' component is not of TENxSpatialList class"
    else
        isSL
    isFL && isSL
}

S4Vectors::setValidity2("TENxVisium", .validTENxVisium)

setMethod("import", "TENxVisium", function(con, format, text, ...) {
    sce <- import(con@resources)
    slist <- import(con@spatialList)
    img <- slist[["imgData"]]
    spd <- slist[["colData"]]
    matches <- intersect(
        colnames(sce),
        rownames(spd)
    )
    spd <- spd[matches, ]
    sce <- sce[, matches]

    SpatialExperiment::SpatialExperiment(
        assays = assays(sce),
        rowData = DataFrame(Symbol = rowData(sce)[["Symbol"]]),
        sample_id = con@spatialList@images,
        colData = spd,
        spatialCoordsNames = con@coordNames,
        imgData = img
    )
})
