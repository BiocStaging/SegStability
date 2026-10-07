#' Descriptive segmentation sensitivity in spatial transcriptomics
#'
#' Infrastructure for comparing segmentations of an immutable molecule universe.
#' Create universe-bound assignments with [SegmentationPartition()] and inspect
#' the full census, semantic statuses and retained ambiguity through its
#'   accessors.
#' Bind named partitions to one shared universe with [SegmentationSet()] and
#' inspect molecule/feature data and segmentation information.
#' Use [compareSegmentations()] for complete shared-molecule correspondence,
#' independent assigned denominators, directional dominance and components
#' including isolates. Inspect results with [correspondenceEdges()],
#' [cellSummary()], [componentSummary()] and [assignmentSummary()].
#' Inspect scientific consistency and integrity with [validateSegmentation()]
#' and [validateSegmentationSet()], metadata with [validateProvenance()], and
#' spatial inputs with [validateCoordinates()] and [validateGeometry()].
#' Validation reports preserve unknown evidence and never repair input.
#'
#' @keywords internal
#' @md
#' @importFrom methods is show
#' @importFrom S4Vectors DataFrame
#' @importFrom digest digest
"_PACKAGE"
