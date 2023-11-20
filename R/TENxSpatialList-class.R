.TENxSpatialList <- setClass(
    "TENxSpatialList",
    contains = "TENxFileList",
    slots = c(
        images = "character",
        scaleJSON = "character",
        tissuePos = "character"
    )
)

.check_file <- function(obj, pattern) {
    fname <- switch(
        pattern,
        "positions.*\\.csv$" = "tissue positions",
        "scalefactors.*\\.json$" = "scalefactor JSON"
    )
    if (any(grepl(pattern, names(obj))))
        TRUE
    else
        paste0("The '", fname, "' file was not found")
}

.validTENxSpatialList <- function(object) {
    .check_file(object, "positions.*\\.csv$")
    .check_file(object, "scalefactors.*\\.json$")
}

S4Vectors::setValidity2("TENxSpatialList", .validTENxSpatialList)

.SCALE_JSON_FILE <- "scalefactors_json.json"

TENxSpatialList <- function(
    resource, images = c("lowres", "hires", "detected", "aligned"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions.*\\.csv"
) {
    images <- match.arg(images, several.ok = TRUE)
    compressed <- endsWith(resource, ".tar.gz")
    spatf <- TENxFileList(resource, compressed = compressed)
    if (compressed)
        spatf <- decompress(con = spatf)
    tissuePos <- grep(tissuePattern, names(spatf), value = TRUE)
    if (!length(tissuePos))
        stop("No tissue positions file found with pattern: ", tissuePattern)

    .TENxSpatialList(
        spatf, images = images, scaleJSON = jsonFile, tissuePos = tissuePos
    )
}

setMethod("import", "TENxSpatialList", function(con, format, text, ...) {
    jsonFile <- con@scaleJSON
    sfs <- jsonlite::fromJSON(txt = path(con[jsonFile]))

    DFs <- lapply(con@images, function(image) {
        .getImgRow(con = con, image = image, scaleFx = sfs)
    })
    list(
        imgData = as(do.call(rbind, DFs), "DataFrame"),
        colData = import(
            TENxSpatialCSV(
                path(con)[con@tissuePos]
            )
        )
    )
})

.getImgRow <- function(con, image, scaleFx) {
    scfactor <- NA_integer_
    fileNames <- names(con)
    filePaths <- path(con)
    imgFile <- grep(image, fileNames, value = TRUE)
    imgPath <- filePaths[endsWith(filePaths, imgFile)]
    if (!length(imgPath))
        stop("Image: ", image, " not found in list of file names")
    spi <- SpatialExperiment::SpatialImage(imgPath)
    scaleName <- grep(image, names(scaleFx), value = TRUE)
    if (length(scaleName))
        scfactor <- unlist(scaleFx[scaleName])

    DataFrame(
        image_id = image,
        data = I(list(spi)),
        scaleFactor = scfactor
    )
}
