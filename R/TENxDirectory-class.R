.TENxFileList <- setClass(
    Class = "TENxFileList",
    contains = "SimpleList"
)

.validTENxFileList <- function(object) {
    validClasses <- vapply(object, function(ifile) {
        is(ifile, "TENxFile")
    }, logical(1L))
    if (all(validClasses))
        TRUE
    else
        "Some files are not of class 'TENxFile'"
}

S4Vectors::setValidity2("TENxFileList", .validTENxFileList)
