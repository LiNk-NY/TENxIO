library(DropletTestFiles)
# listTestFiles()
library(ExperimentHub)
eh <- ExperimentHub()

path <- eh[["EH3686"]]
expect_warning(
    res <- TENxH5(
        resource = path, group = "mm10", ranges = NA_character_,
        version = "2", extension = "h5", rowidx = 1:10, colidx = 1:10
    ),
    pattern = "'group' not in known.*"
)
expect_true(
    is(res, "TENxH5")
)
expect_true(
    is(
        import(res, ref = NA),
        "SingleCellExperiment"
    )
)
expect_true(
    is(
        import(res),
        "SingleCellExperiment"
    )
)

path <- eh[["EH3688"]]
expect_error(
    TENxFileList(path, compressed = TRUE),
    "No extension present"
)
expect_true(
    is(
        res <- TENxFileList(
            path, extension = "tar.gz", compressed = TRUE, version = "2"
        ),
        "TENxFileList"
    )
)
expect_true(
    is(
        import(res),
        "SingleCellExperiment"
    )
)
