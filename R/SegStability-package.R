#' Descriptive segmentation sensitivity in spatial transcriptomics
#'
#' Infrastructure for comparing segmentations of an immutable molecule universe.
#' Create universe-bound assignments with [SegmentationPartition()] and inspect
#' the full census, semantic statuses and retained ambiguity through its
#'   accessors.
#' Bind named partitions to one shared universe with [SegmentationSet()] and
#' inspect molecule/feature data and segmentation information.
#' Comparative scientific analysis functions are not implemented at this stage.
#'
#' @keywords internal
#' @importFrom methods is show
#' @importFrom S4Vectors DataFrame
#' @importFrom digest digest
"_PACKAGE"
