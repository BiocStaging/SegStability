#' Descriptive segmentation sensitivity in spatial transcriptomics
#'
#' SegStability compares alternative cell segmentations of one immutable
#' molecule universe. [SegmentationPartition()] records each full cell census
#' and the assigned, unassigned or ambiguous status of every molecule.
#' [SegmentationSet()] binds named partitions to the shared universe; retrieve
#' one with `ss[["name"]]` and read keyed assignments with [assignmentData()].
#'
#' [compareSegmentations()] builds the complete shared-molecule graph with
#' independent assigned denominators. Its edge, cell, component and status
#' summaries retain empty cells and exact dominance ties. [stateStability()]
#' uses separately supplied categorical labels. [supportSensitivity()]
#' describes directional edge retention without changing the primary graph.
#' Optional [anchorTopology()] and [compareGeometry()] use qualified spatial
#' evidence; [replicateAgreement()] compares runs through molecule identity.
#'
#' [fromMoleculeExperiment()] and
#' [geometryFromSpatialFeatureExperiment()] accept explicit IDs and mappings
#' from optional Bioconductor containers. Validation functions report missing
#' evidence without inferring provenance or repairing inputs. Results are
#' descriptive and do not rank segmentation accuracy or biological truth.
#'
#' @examples
#' data("segmentationExample", package = "SegStability")
#' names(segmentationExample)
#' head(segmentationExample$assignments$first, 3)
#' @seealso [segmentationExample]
#' @md
#' @importFrom methods is show
#' @importFrom S4Vectors DataFrame
#' @importFrom digest digest
"_PACKAGE"
