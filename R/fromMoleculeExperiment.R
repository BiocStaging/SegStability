#' Import retained molecules with explicit segmentation evidence
#'
#' Translate one MoleculeExperiment assay/sample into constructor arguments
#' for descriptive segmentation sensitivity. Original identity must survive
#' upstream processing; coordinates or row order cannot recover lost identity.
#'
#' @param object A MoleculeExperiment with retained original molecule IDs.
#' @param assay Explicit character assay name.
#' @param sample_id Explicit character source sample key; not an inferred
#'   specimen identifier.
#' @param name Segmentation name.
#' @param specimen_id Explicit specimen/ROI identity.
#' @param cells Complete data.frame or DataFrame cell census with cell_id,
#'   including zero-molecule cells.
#' @param assignments data.frame or DataFrame containing molecule_id, status
#'   and cell_id. Status is explicitly assigned, unassigned or ambiguous.
#' @param molecule_id_col Name of the retained original molecule ID column.
#' @param coordinate_cols Named character vector mapping x/y or x/y/z to source
#'   columns, or NULL to disable coordinate translation.
#' @param universe Optional canonical molecule universe data.frame/DataFrame.
#'   Exact typed membership, feature identity and shared columns must agree.
#' @param ambiguity External-ID ambiguity evidence as documented by
#'   [SegmentationPartition()]. NULL retains unknown candidate evidence.
#' @param coordinate_info Coordinate declarations as a SimpleList. Unknown
#'   registration stays unknown; coordinates are never transformed.
#' @param run_metadata Explicit run metadata SimpleList.
#' @param provenance Explicit provenance SimpleList; no values are inferred.
#' @param metadata Ancillary metadata SimpleList, retained in partition and set.
#' @return A SimpleList of exactly the six [SegmentationSet()] arguments,
#'   containing a single validated partition. Use
#'   `do.call(SegmentationSet, as.list(result))` to construct the set.
#' @details Feature identity comes from retained character feature keys.
#'   Extraction follows table order; IDs never come from that order. Explicit
#'   assignments and a full census are required because upstream boundaries
#'   and assigned molecule rows cannot establish completeness or ambiguity.
#'   Default upstream essential imports may discard molecule IDs and are then
#'   unsupported. All certified ID families are preserved, including optional
#'   integer64. Source tables require compatible schemas; extra columns survive.
#'   Canonical ID/feature/coordinate columns precede ancillary columns. With a
#'   supplied universe, its column order is retained instead.
#'   Source object, geometry and assays are not modified or recomputed.
#' @examples
#' if (requireNamespace("MoleculeExperiment", quietly = TRUE)) {
#'     tab <- data.frame(molecule_id = c("001", "002"),
#'         x_location = c(1, 2), y_location = c(3, 4))
#'     me <- MoleculeExperiment::MoleculeExperiment(
#'         molecules = list(detected = list(roi = list(g = tab))))
#'     args <- fromMoleculeExperiment(me, "detected", "roi", "native", "roi1",
#'         cells = data.frame(cell_id = c("cell1", "empty")),
#'         assignments = data.frame(molecule_id = c("002", "001"),
#'             status = c("unassigned", "assigned"),
#'             cell_id = c(NA_character_, "cell1")))
#'     ss <- do.call(SegmentationSet, as.list(args))
#'     moleculeData(ss)
#' }
#' @export
fromMoleculeExperiment <- function(object, assay, sample_id, name, specimen_id,
    cells, assignments, molecule_id_col = "molecule_id",
    coordinate_cols = c(x = "x_location", y = "y_location"), universe = NULL,
    ambiguity = NULL, coordinate_info = S4Vectors::SimpleList(),
    run_metadata = S4Vectors::SimpleList(),
    provenance = S4Vectors::SimpleList(), metadata = S4Vectors::SimpleList()) {
    .adapterPackage("MoleculeExperiment", object)
    .adapterScalar(assay, "assay")
    .adapterScalar(sample_id, "sample_id")
    .adapterScalar(name, "name")
    .adapterScalar(specimen_id, "specimen_id")
    .adapterScalar(molecule_id_col, "molecule_id_col")
    .adapterSpecimen(coordinate_info, specimen_id)
    molecules <- .meMolecules(object, assay, sample_id, molecule_id_col,
                              coordinate_cols)
    molecules <- .meUniverse(molecules, universe)
    .adapterSourceSpecimen(molecules, specimen_id)
    .adapterValidation(validateCoordinates(molecules, coordinate_info))
    cells <- .ingestionTable(cells, "cells")
    assignments <- .ingestionTable(assignments, "assignments")
    if (anyDuplicated(names(assignments)) ||
        !all(c("molecule_id", "status", "cell_id") %in% names(assignments))) {
        stop("assignments require molecule_id, status and cell_id",
             call. = FALSE)
    }
    partition <- SegmentationPartition(molecules, name, cells,
        assignments$molecule_id, assignments$status, assignments$cell_id,
        ambiguity, run_metadata, provenance, metadata)
    segmentations <- S4Vectors::SimpleList(partition)
    names(segmentations) <- name
    S4Vectors::SimpleList(molecules = molecules, segmentations = segmentations,
        specimen_id = specimen_id, coordinate_info = coordinate_info,
        provenance = provenance, metadata = metadata)
}
