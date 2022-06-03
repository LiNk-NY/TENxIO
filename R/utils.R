.setNames <- function(object = nm, nm) {
    names(object) <- nm
    object
}

.checkPkgsAvail <- function(pkgnames, walkback = -3L) {
    vapply(pkgnames, function(pkgname) {
        func <- as.character(sys.call(walkback)[[1L]])
        func <- tail(func, 1L)
        if (!requireNamespace(pkgname, quietly = TRUE))
            stop("Install '", pkgname, "' to use '", func, "'", call. = FALSE)
        else
            TRUE
    }, logical(1L))
}
