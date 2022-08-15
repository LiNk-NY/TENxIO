fi <- system.file(
    "extdata", "pbmc_granulocyte_sorted_3k_ex_atac_peak_annotation.tsv",
    package = "TENxIO", mustWork = TRUE
)

peak_file <- TENxPeaks(fi)

expect_true(
    is(peak_file, "TENxPeaks")
)

expect_identical(
    peak_file@extension, "tsv"
)

peak_anno <- import(peak_file)

expect_true(
    is(peak_anno, "GRanges")
)

expect_identical(
    length(peak_anno), 10L
)

expect_identical(
    colnames(mcols(peak_anno)), c("gene", "distance", "peak_type")
)
