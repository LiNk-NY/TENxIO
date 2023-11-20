.TENxSpatialList <- setClass(
    "TENxSpatialList",
    contains = "TENxFileList",
    slots = c(
        image = "character"
    )
)

.check_file <- function(obj, pattern) {
    if (any(grepl(pattern, names(obj))))
        TRUE
    else
        paste0("The '", name, "' file was not found")
}

.validTENxSpatialList <- function(object) {
    .check_file(object, "positions.*\\.csv$")
    .check_file(object, "scalefactors.*\\.json$")
}

S4Vectors::setValidity2("TENxSpatialList", .validTENxSpatialList)

TENxSpatialList <- function(
    resource, image = c("lowres", "hires", "detected", "aligned")
) {
    image <- match.arg(image, several.ok = TRUE)
    compressed <- endsWith(resource, ".tar.gz")
    spat <- TENxFileList(resource, compressed = compressed)
    spatf <- decompress(con = spat)
    imgs <- grep(paste(image, collapse = "|"), names(spatf), value = TRUE)
    
    .TENxSpatialList(
        spatf, image = imgs
    )
}