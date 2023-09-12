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
