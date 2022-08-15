h5f <- system.file(
    "extdata", "pbmc_granulocyte_ff_bc_ex.h5",
    package = "TENxIO", mustWork = TRUE
)
tenxh5 <- TENxH5(h5f)
expect_true(
    is(tenxh5, "TENxH5")
)
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
