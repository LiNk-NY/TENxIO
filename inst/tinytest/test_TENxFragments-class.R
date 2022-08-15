fr <- system.file(
    "extdata", "pbmc_3k_atac_ex_fragments.tsv.gz",
    package = "TENxIO", mustWork = TRUE
)
## yieldSize warning
expect_warning(
    tfr <- TENxFragments(fr)
)

tfr <- TENxFragments(fr, extension = "tsv.gz")
expect_identical(
    tfr@extension, "tsv.gz"
)
expect_silent(
    tfr <- TENxFragments(fr, yieldSize = 250)
)

expect_true(
    is(tfr, "TENxFragments")
)

fra <- import(tfr)
expect_true(
    is(fra, "RaggedExperiment")
)

expect_identical(
    assayNames(fra), c("barcode", "readSupport")
)

expect_identical(
    dim(fra), c(10L, 10L)
)
