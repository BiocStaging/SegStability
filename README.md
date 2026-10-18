# SegStability

SegStability compares alternative cell segmentations of the same
molecule-resolved spatial transcriptomics experiment. It counts how exact
molecule identities move between cells and describes how those changes affect
categorical cell-state conclusions.

## Why use it?

Two boundary sets can assign different cells to a molecule, leave it
unassigned, or mark it ambiguous. SegStability keeps those outcomes separate.
Its weighted correspondence graph uses a complete cell census for each
segmentation, including empty cells, and independent assigned-molecule
denominators.

## Installation

Version 0.99.0 is a pre-submission development version. It is not yet
available through Bioconductor's package repository. With R and the matching
Bioconductor devel release installed, the public source repository can be
installed with:

```r
BiocManager::install(c("S4Vectors", "BiocStyle"))
remotes::install_github("ramakrishna-p21/SegStability", build_vignettes = TRUE)
```

The `remotes` package is needed for `install_github()`; install it first if
necessary. Optional spatial and container adapters have their own suggested
dependencies.

## A small example

The bundled tables describe 24 synthetic molecules and two segmentations.
No biological interpretation is attached to their coordinates or labels.

```r
library(SegStability)
data("segmentationExample", package = "SegStability")
d <- segmentationExample

make_partition <- function(side) {
    a <- d$assignments[[side]]
    ambiguity <- d$ambiguity[[side]]
    if (!is.null(ambiguity)) {
        ambiguity <- S4Vectors::SimpleList(
            evidence = S4Vectors::DataFrame(ambiguity$evidence),
            candidates = S4Vectors::DataFrame(ambiguity$candidates))
    }
    SegmentationPartition(d$molecules, side, d$cells[[side]],
                          a$molecule_id, a$status, a$cell_id,
                          ambiguity = ambiguity)
}

ss <- SegmentationSet(d$molecules,
    list(first = make_partition("first"),
         second = make_partition("second")),
    specimen_id = "synthetic-example")
cmp <- compareSegmentations(ss, "first", "second")
componentSummary(cmp)[, c("component_id", "component_class")]
stateStability(cmp, d$states)$global
```

Retrieve a stored partition with `ss[["first"]]`, its keyed assignments with
`assignmentData(ss, "first")`, and stored coordinate declarations with
`summary(ss)$coordinate_info` or `summary(cmp)$coordinate_info`.

## Further analyses

`supportSensitivity()` reports which original shared-molecule edges remain at
directional support thresholds, while keeping the original denominators.
Optional `anchorTopology()` and `compareGeometry()` use qualified spatial
inputs; `replicateAgreement()` compares stochastic segmentation runs on a
shared molecule universe. MoleculeExperiment and SpatialFeatureExperiment
adapters require explicit identity, census, and geometry mappings.

These results describe sensitivity. They do not estimate segmentation
accuracy, recover lost molecule IDs, or select a biologically correct
segmentation.

## Documentation and support

The [vignette](vignettes/SegStability.Rmd) develops the workflow and explains
the output tables. Installed help is available through `?SegStability` and
`?compareSegmentations`. Report bugs and enhancement requests in
[GitHub Issues](https://github.com/ramakrishna-p21/SegStability/issues).
For general usage questions, use the
[Bioconductor Support Site](https://support.bioconductor.org/) and include
`SegStability` as a tag once the package is listed there.

For citation details, run `citation("SegStability")` after installation.
