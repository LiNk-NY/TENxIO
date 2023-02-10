## Changes in version 1.0.1

### Bug fixes and minor improvements

* The `import` method for `TENxFileList` was returning a nested
`SummarizedExperiment` within the `SingleCellExperiment`. The `counts` and
`assay(..., i="counts")` methods should only return the bare `Matrix` rather
than the embedded `SummarizedExperiment`.

## Changes in version 1.0.0

* Package released in Bioconductor!
