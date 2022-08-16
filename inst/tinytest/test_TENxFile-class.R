expect_error(
    TENxFile(tempfile())
)

test_file <- tempfile()
file.create(test_file)

tenxfile <- TENxFile(test_file)
expect_true(
    is(tenxfile, "TENxFile")
)

expect_identical(
    tenxfile@resource, test_file
)

expect_identical(
    tenxfile@extension, ""
)

file.remove(test_file)

test_file <- tempfile(fileext = ".bed")
file.create(test_file)

tenxfile <- TENxFile(test_file)

expect_identical(
    tenxfile@extension, "bed"
)

pbmc_url <-
    "https://raw.githubusercontent.com/waldronlab/TENxIO/devel/inst/extdata/10k_pbmc_ATACv2_f_bc_ex.h5"

remoteh5 <- TENxFile(pbmc_url)

expect_true(
    is(remoteh5, "TENxH5")
)

expect_true(
    remoteh5@remote
)
