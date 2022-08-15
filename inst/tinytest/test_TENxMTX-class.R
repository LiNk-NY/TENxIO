mtxf <- system.file(
    "extdata", "pbmc_3k_ff_bc_ex.mtx",
    package = "TENxIO", mustWork = TRUE
)

con <- TENxMTX(mtxf)

expect_true(
    is(con, "TENxMTX")
)

expect_identical(
    con@extension, "mtx"
)

expect_false(
    con@compressed
)

mtx <- import(con)
expect_true(
    is(mtx, "SummarizedExperiment")
)

expect_identical(
    assayNames(mtx), "counts"
)

expect_identical(
    dim(mtx), c(171L, 10L)
)
