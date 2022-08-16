fl <- system.file(
    "extdata", "pbmc_granulocyte_sorted_3k_ff_bc_ex_matrix.tar.gz",
    package = "TENxIO", mustWork = TRUE
)

txfl <- TENxFileList(fl)

expect_identical(
    length(txfl), 1L
)

expect_true(
    txfl@compressed
)

expect_identical(
    txfl@extension, "tar.gz"
)

fli <- import(txfl)

expect_true(
    is(fli, "SingleCellExperiment")
)

expect_identical(
    dim(fli), c(10L, 10L)
)

expect_identical(
    names(rowData(fli)), c("ID", "Symbol", "Type")
)

expect_identical(
    mainExpName(fli), "Gene Expression"
)

expect_identical(
    assayNames(fli), "counts"
)
