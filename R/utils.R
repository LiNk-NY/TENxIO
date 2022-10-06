.setNames <- function(object = nm, nm) {
    names(object) <- nm
    object
}

.checkPkgsAvail <- function(pkgnames) {
    toinst <- pkgnames[!pkgnames %in% rownames(utils::installed.packages())]
    if (length(toinst))
        .install_suggestion(toinst)
    else
        TRUE
}

.install_suggestion <- function(pkgs) {
    n <- length(pkgs)
    txt <- if (identical(n, 1L)) '"%s"' else 'c(\n    "%s"\n  )'
    fmt <-'  BiocManager::install(%s)'
    fmt <- sprintf(fmt, txt)
    pkgs <- paste(strwrap(
        paste(pkgs, collapse='", "'),
        width = getOption("width") - 4L
    ), collapse="\n    ")
    stop("Enable this functionality with",
        "\n\n", sprintf(fmt, pkgs), "\n\n",
        sep = "", call. = FALSE
    )
}
