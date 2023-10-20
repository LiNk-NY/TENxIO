#' @docType class
#'
#' @title A class to represent and import Visium data
#'
#' @description This class is a composed class of [TENxFileList] which can
#'   contain a list of [TENxFile] objects and a [TENxSpatialList] object. It is
#'   meant to handle Visium data from 10X Genomics.
#'
#' @details Typically, the user will not create an object of this class directly
#'   but rather use the [TENxVisium()] constructor function to create an object
#'   of this class.
#'
#' @slot resources A [TENxFileList] object containing the Visium data.
#'
#' @slot spatialList A [TENxSpatialList] object containing the spatial
#'
#' @slot coordNames `character()` A vector specifying the names
#'   of the columns in the spatial data containing the spatial coordinates.
#'
#' @slot sampleId `character(1)` A scalar specifying the sample identifier.
#'
#' @return A [SpatialExperiment] object
#'
#' @exportClass TENxVisium
.TENxVisium <- setClass(
    Class = "TENxVisium",
    slots = c(
        resources = "TENxFileList",
        spatialList = "TENxSpatialList",
        coordNames = "character",
        sampleId = "character"
    )
)

#' @rdname TENxVisium-class
#'
#' @title Represent and import Visium data from 10X Genomics
#'
#' @description `TENxVisium` is a class to represent and import Visium data. It
#'   is a composed class of [TENxFileList] which can contain a list of
#'   [TENxFile] objects and a [TENxSpatialList] object.
#'
#' @param resources A [TENxFileList] object or a file path to the tarball
#'   containing the matrix / assay data resources.
#'
#' @param spatialResource A [TENxSpatialList] object or a file path to the
#'   tarball containing the spatial data.
#'
#' @param sample_id `character(1)` A single string specifying the sample ID.
#'
#' @param images `character()` A vector specifying the images to be imported;
#'   can be one or multiple of "lowres", "hires", "detected", "aligned".
#'
#' @param jsonFile `character(1)` A single string specifying the name of the
#'  JSON file containing the scale factors.
#'
#' @param tissuePattern `character(1)` A single string specifying the pattern
#'   to match the tissue positions file.
#'
#' @param spatialCoordsNames `character()` A vector of strings specifying the
#'  names of the columns in the spatial data containing the spatial coordinates.
#'
#' @param ... In the constructor, additional arguments passed to
#'   [TENxFileList()]; otherwise, not used.
#'
#' @examples
#' \dontrun{
#'     spatialtar <- "~/data/V1_Adult_Mouse_Brain_spatial.tar.gz" 
#'     dir.create(sdir <- tempfile())
#'     untar(spatialtar, exdir = sdir)
#'
#'     matrixtar <-
#'         "~/data/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz" 
#'     dir.create(mdir <- tempfile())
#'     untar(matrixtar, exdir = mdir)
#'
#'     tv <- TENxVisium(
#'         resources = mdir, spatialResource = tdir, images = "lowres"
#'     )
#'     import(tv)
#'
#'     unlink(sdir, recursive = TRUE)
#'     unlink(mdir, recursive = TRUE)
#' }
#' @export
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

#' @describeIn TENxVisium-class Import Visium data
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom BiocIO import
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SummarizedExperiment assays rowData colData
#'
#' @exportMethod import
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
        rowData = S4Vectors::DataFrame(Symbol = rowData(sce)[["Symbol"]]),
        sample_id = con@sampleId,
        colData = spd,
        spatialCoordsNames = con@coordNames,
        imgData = img
    )
})
