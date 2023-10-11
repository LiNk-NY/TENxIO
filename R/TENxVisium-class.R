.TENxVisium <- setClass(
    Class = "TENxVisium",
    slots = c(
        resources = "TENxFileList",
        spatialList = "TENxSpatialList",
        coordNames = "character",
        sampleId = "character"
    )
)

TENxVisium <- function(
    resources,
    spatialResource,
    sample_id = "sample01",
    images = c("lowres", "hires", "detected", "aligned"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions.*\\.csv",
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    ...
) {
    images <- match.arg(images, several.ok = TRUE)
    if (!is(resources, "TENxFileList"))
        resources <- TENxFileList(resources, ...)
    if (!is(spatialResource, "TENxSpatialList"))
        spatialResource <- TENxSpatialList(
            resource = spatialResource, sample_id = sample_id, images = images,
            jsonFile = jsonFile, tissuePattern = tissuePattern
        )

    .TENxVisium(
        resources = resources,
        spatialList = spatialResource,
        coordNames = spatialCoordsNames,
        sampleId = sample_id
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
        sample_id = con@sampleId,
        colData = spd,
        spatialCoordsNames = con@coordNames,
        imgData = img
    )
})
