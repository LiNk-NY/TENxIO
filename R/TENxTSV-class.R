#' TENxTSV: A class to represent 10x tab separated values files
#'
#' @aliases TENxTSV
#'
#' @description This class is general purpose for reading in tabular data from
#'   the 10x Genomics website with the `.tsv` file extension. The class also
#'   supports compressed files, i.e., those with the `.tsv.gz` extension.
#'
#' @details Typical `.tsv` files obtained from the 10X website are compressed
#'   and contain information relevant to 'barcodes' and 'features'. Currently,
#'   the code only supports files such as `features.tsv.*` and `barcodes.tsv.*`.
#'
#' @inheritParams TENxMTX-class
#'
#' @return A `TENxTSV` class object; a `tibble` for the import method
#'
#' @exportClass TENxTSV
.TENxTSV <- setClass(Class = "TENxTSV", contains = "TENxFile")

#' @describeIn TENxTSV-class General import method for `tsv` files from 10x;
#'   using `readr::read_tsv` and returning a `tibble` representation
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom readr read_tsv
#'
#' @export
setMethod("import", "TENxTSV", function(con, format, text, ...) {
    resource <- path(con)
    if (identical(Sys.info()[["sysname"]], "Darwin"))
        readr::local_edition(1)
    df <- readr::read_tsv(
        resource, col_names = FALSE, show_col_types = FALSE, progress = FALSE,
        ...
    )
    fname <- basename(resource)
    if (startsWith(fname, "features.tsv"))
        names(df) <- c("ID", "Symbol", "Type", "Chr", "Start", "End")
    else if (startsWith(fname, "barcodes.tsv"))
        names(df) <- "barcode"
    df
})

#' @rdname TENxTSV-class
#'
#' @inheritParams TENxMTX
#'
#' @export
TENxTSV <-  function(resource, compressed = FALSE, ...) {
    dots <- list(...)
    ext <- dots[["extension"]]
    if (is.null(ext))
        ext <- .get_ext(resource)
    compr <- identical(ext, "tsv.gz")
    if (!ext %in% c("tsv.gz", "tsv"))
        warning("File extension is not 'tsv'; import may fail", call. = FALSE)
    .TENxTSV(resource = resource, compressed = compr, extension = ext)
}
