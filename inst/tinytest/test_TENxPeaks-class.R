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

sce <- SingleCellExperiment::SingleCellExperiment(
    assays = SimpleList(counts = matrix(rnorm(100), 10, 10))
)

expect_error(
    annotation(sce, name = "peak_annotation") <- peak_anno
)

expect_silent(
    annotation(sce, name = "peak_annotation") <- peak_file
)

expect_identical(
    metadata(sce), list(annotation = list(peak_annotation = peak_anno))
)

expect_identical(
    names(annotation(sce)), "peak_annotation"
)

expect_warning(
    annotation(sce, append = FALSE) <- peak_file
)

expect_silent(
    annotation(sce, append = TRUE) <- peak_file
)

expect_identical(
    length(metadata(sce)), 1L
)

expect_identical(
    length(metadata(sce)[["annotation"]]), 2L
)

expect_identical(
    names(metadata(sce)[["annotation"]]),
    c("peak_annotation", "peak_annotation")
)

