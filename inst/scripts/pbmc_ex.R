## subset H5 file
library(rhdf5)
# setwd("~/data/10x/pbmc_3k")
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
as.object_size(file.info(h5new)$size, "kB")
#' [1] "30.2 kB"


h5f <- "10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5"
h5ls(h5f)
newmat <- HDF5Array::TENxMatrix(h5f, "matrix")[1:10, 1:10]
h5new <- "inst/extdata/10k_pbmc_ATACv2_f_bc_ex.h5"
# file.remove(h5new)
HDF5Array::writeTENxMatrix(newmat, h5new, group = "matrix", verbose = TRUE)

tkeys <- h5read(h5f, "/matrix/features/_all_tag_keys") ## same
myidx <- list(1:10)
h1 <- h5read(h5f, "/matrix/features/feature_type", index = myidx)
h2 <- h5read(h5f, "/matrix/features/genome", index = myidx)
h3 <- h5read(h5f, "/matrix/features/id", index = myidx)
h4 <- h5read(h5f, "/matrix/features/derivation", index = myidx)
h5 <- h5read(h5f, "/matrix/features/name", index = myidx)

h5createGroup(h5new, "/matrix/features/")
h5write(tkeys, h5new, "/matrix/features/_all_tag_keys")
h5write(h1, h5new, "/matrix/features/feature_type")
h5write(h2, h5new, "/matrix/features/genome")
h5write(h3, h5new, "/matrix/features/id")
h5write(h4, h5new, "/matrix/features/derivation")
h5write(h5, h5new, "/matrix/features/name")

h5ls(h5new)
aa <- TENxH5(h5new, ranges = "/features/id")
import(aa)

as.object_size(file.size(h5new), "kB")
#' [1] "30.1 kB"


R.utils::gunzip(
    "filtered_feature_bc_matrix/matrix.mtx.gz",
    "filtered_feature_bc_matrix/test_matrix.mtx"
)
gg <- readLines("test_matrix.mtx", n = 13)
gg[[3]] <- "171 10 10"

writeLines(gg, "test_write_matrix.mtx")
assay(import(TENxMTX("test_write_matrix.mtx")))

file.copy("test_write_matrix.mtx", "inst/extdata/pbmc_3k_ff_bc_ex.mtx")

as.object_size(file.size(mtxf), unit = "kB")
#' [1] "0.2 kB"


# TENxFragments -----------------------------------------------------------

library(Rsamtools)
fr <- "pbmc_granulocyte_sorted_3k_atac_fragments.tsv"
tf <- TabixFile(fr, yieldSize = 10)
examp <- read.table(textConnection(scanTabix(tf)[[1]]))

names(examp) <- c("chrom", "chromStart", "chromEnd", "barcode", "readSupport")
readr::write_tsv(
    examp, "inst/extdata/pbmc_3k_atac_ex_fragments.tsv", col_names = FALSE
)
#' [1] "0.2 kB"
system2("bgzip inst/extdata/pbmc_3k_atac_ex_fragments.tsv", stdout = TRUE)
system("tabix -p bed inst/extdata/pbmc_3k_atac_ex_fragments.tsv.gz")

fragf <- "inst/extdata/pbmc_3k_atac_ex_fragments.tsv.gz"
as.object_size(file.size(fragf), unit = "kB")
#' [1] "0.2 kB"
as.object_size(file.size(paste0(fragf, ".tbi")), unit = "kB")
#' [1] "0.1 kB"

aa <- import(TENxFragments(fragf, yieldSize = 10))


# Peaks -------------------------------------------------------------------

fi <- "pbmc_granulocyte_sorted_3k_atac_peak_annotation.tsv"
pkt <- readr::read_tsv(fi, n_max = 10)
newpeak <- "inst/extdata/pbmc_granulocyte_sorted_3k_ex_atac_peak_annotation.tsv"
readr::write_tsv(pkt, newpeak)
as.object_size(file.size(newpeak), unit = "kB")
#' [1] "0.5 kB"

import(TENxPeaks(newpeak))


# TENxFileList ------------------------------------------------------------

basedir <- "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix/filtered_feature_bc_matrix/"
bc <- readr::read_tsv(
    file.path(basedir, "barcodes.tsv.gz"), col_names = FALSE, n_max = 10
)
readr::write_tsv(
    bc, file.path(basedir, "barcodes.tsv.gz"), col_names = FALSE
)
ft <- readr::read_tsv(
    file.path(basedir, "features.tsv.gz"), col_names = FALSE, n_max = 10
)
readr::write_tsv(
    ft, file.path(basedir, "features.tsv.gz"), col_names = FALSE
)
# file.copy("TENxIO/inst/extdata/pbmc_3k_ff_bc_ex.mtx",
#   file.path(basedir, "matrix.mtx"))
R.utils::gzip(
    file.path(basedir, "matrix.mtx"), file.path(basedir, "matrix.mtx.gz"),
    overwrite = TRUE
)
tar(
    "TENxIO/inst/extdata/pbmc_granulocyte_sorted_3k_ff_bc_ex_matrix.tar.gz",
    "./pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix",
    compression = "gzip"
)

import(
    TENxFileList(
    "TENxIO/inst/extdata/pbmc_granulocyte_sorted_3k_ff_bc_ex_matrix.tar.gz"
    )
)

as.object_size(
    file.size("pbmc_granulocyte_sorted_3k_ff_bc_ex_matrix.tar.gz"), "KB"
)
#' [1] "0.7 kB"
