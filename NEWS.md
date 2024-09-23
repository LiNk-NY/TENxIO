## Changes in version 1.8.0

### New features

* The `import` method for h5 files now includes all `rowData` which typically
includes "ID", "Symbol", and "Type" columns.

### Bug fixes and minor improvements

* `rownames` are set to first column in the `features.tsv.gz` data (`rownames`
for h5 files are determined by `HDF5Array::TENxMatrix`)
* `TENxH5` now tests whether there is a `/matrix/features/interval` dataset in
the h5 file. It sets ranges to `NA_character_` when interval data is not found.

## Changes in version 1.6.0

### Bug fixes and minor improvements

* Use `rowData` from filtered `SingleCellExperiment` object (@MalteThodberg, #5)
* Add examples to `TENxFileList` documentation and improve constructor function

## Changes in version 1.4.0

### Bug fixes and minor improvements

* Skip unit tests when remote H5 access is not configured.

## Changes in version 1.2.0

### New features

* `TENxTSV` class has been added to handle compressed and uncompressed TSV
files.

### Bug fixes and minor improvements

* The `import` method for `TENxFileList` was returning a nested
`SummarizedExperiment` within the `SingleCellExperiment`. The `counts` and
`assay(..., i="counts")` methods should only return the bare `Matrix` rather
than the embedded `SummarizedExperiment` (@dgastn, #2).

## Changes in version 1.0.0

* Package released in Bioconductor!
