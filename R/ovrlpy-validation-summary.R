#' Public summary of the qualified ovrlpy/Xenium application
#'
#' The small CSV at `inst/extdata/validation/ovrlpy_validation_summary.csv`
#' records source identifiers, exact file hashes, the source-defined molecule
#' universe and ROI, and checked vendor/Baysor correspondence and directional
#' support counts. Its `key`, `value`, and `evidence` columns distinguish facts
#' checked in the public-source replay from retained qualification evidence.
#'
#' The evaluated region is tile (9,9) of replicate 1 of
#' `Xenium_V1_FF_Mouse_Brain_MultiSection_1`. Original source molecule IDs and
#' a complete, independent ROI cell census were preserved on each side.
#' Native cells can cross the ROI boundary; their reported denominators here
#' count their own assigned ROI molecules. Vendor and saved Baysor v0.7.1
#' outputs differed in upstream preparation. The record describes
#' correspondence and support, not an isolated boundary effect, accuracy
#' ranking, or biological ground truth.
#'
#' The optional development-time replay uses
#' `inst/scripts/prepare_ovrlpy_validation.py` followed by
#' `inst/scripts/reproduce_ovrlpy_validation.R`. They require explicitly
#' acquired, hash-checked public source files and run only when invoked by a
#' user. No source molecules, large saved segmentations, or benchmark work
#' files are bundled. The source dataset and Zenodo deposit are credited to
#' their respective creators under their declared CC BY 4.0 terms.
#'
#' @format A CSV with one record per key and columns `key`, `value`, and
#'   `evidence`. Read values as text to preserve identifiers and hashes.
#' @return This documentation topic returns no object. Read the installed CSV
#'   with `system.file()` and `read.csv()`, as shown in the example.
#' @source Tiesmeyer et al., *Identifying 3D signal overlaps in spatial
#'   transcriptomics data with ovrlpy*, Nature Biotechnology (2026),
#'   doi:10.1038/s41587-026-03004-8. Saved Baysor data:
#'   doi:10.5281/zenodo.16736589 (record version 2). Original Xenium output:
#'   `Xenium_V1_FF_Mouse_Brain_MultiSection_1`, 10x Genomics Xenium 1.0.2.
#' @examples
#' path <- system.file("extdata", "validation",
#'     "ovrlpy_validation_summary.csv", package = "SegStability")
#' record <- read.csv(path, colClasses = "character")
#' record[record$key %in% c("roi_universe_size", "shared_mass"), ]
#' @name ovrlpy_validation_summary
#' @md
NULL
