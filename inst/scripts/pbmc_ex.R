## subset H5 file
library(rhdf5)
h5f <- "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
as.object_size(file.info(h5f)$size)
#' [1] "38.8 MB"

h5new <- "inst/extdata/pbmc_granulocyte_ff_bc_ex.h5"
newmat <- HDF5Array::TENxMatrix(h5f, "matrix")[1:10, 1:10]
HDF5Array::writeTENxMatrix(newmat, h5new, group = "matrix", verbose = TRUE)

tkeys <- h5read(h5f, "/matrix/features/_all_tag_keys") ## same
myidx <- list(1:10)
h1 <- h5read(h5f, "/matrix/features/feature_type", index = myidx)
h2 <- h5read(h5f, "/matrix/features/genome", index = myidx)
h3 <- h5read(h5f, "/matrix/features/id", index = myidx)
h4 <- h5read(h5f, "/matrix/features/interval", index = myidx)
h5 <- h5read(h5f, "/matrix/features/name", index = myidx)

h5createGroup(h5new, "/matrix/features/")
h5write(tkeys, h5new, "/matrix/features/_all_tag_keys")
h5write(h1, h5new, "/matrix/features/feature_type")
h5write(h2, h5new, "/matrix/features/genome")
h5write(h3, h5new, "/matrix/features/id")
h5write(h4, h5new, "/matrix/features/interval")
h5write(h5, h5new, "/matrix/features/name")

TENxH5(h5new)
as.object_size(file.info(h5new)$size, "MB")
#' [1] "30.2 kB"
