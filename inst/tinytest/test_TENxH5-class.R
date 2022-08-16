h5f <- system.file(
    "extdata", "pbmc_granulocyte_ff_bc_ex.h5",
    package = "TENxIO", mustWork = TRUE
)
tenxh5 <- TENxH5(h5f)
expect_true(
    is(tenxh5, "TENxH5")
)

expect_identical(
    genome(tenxh5),
    c(chr1 = "GRCh38")
)

# set to NA_character_ to test
tenxh5@ranges <- NA_character_
expect_error(
    genome(tenxh5)
)
expect_error(
    rowRanges(tenxh5)
)
tenxh5@ranges <- "/features/interval"

sceh5 <- import(tenxh5)
expect_true(
    is(sceh5, "SingleCellExperiment")
)

expect_identical(
    mainExpName(sceh5), "Gene Expression"
)

expect_identical(
    colnames(rowData(sceh5)), c("ID", "Symbol", "Type")
)

expect_identical(
    assayNames(sceh5), "counts"
)

expect_identical(
    dim(sceh5), c(10L, 10L)
)
