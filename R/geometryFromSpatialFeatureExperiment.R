#' Extract an explicitly keyed SpatialFeatureExperiment geometry sidecar
#'
#' Preserve existing polygons and CRS for optional spatial analysis. Geometry
#' supplies no evidence of molecule assignment identity or registration.
#'
#' @param object A SpatialFeatureExperiment.
#' @param type Explicit character geometry name.
#' @param sample_id Explicit character source sample key.
#' @param geometry_kind Either "col" or "annot", without automatic fallback.
#' @param cell_map data.frame/DataFrame with unique character geometry_id source
#'   keys and lossless cell_id values. Every selected geometry must be mapped.
#' @param cells Complete census data.frame/DataFrame or SegmentationPartition.
#' @param coordinate_info Explicit coordinate declarations SimpleList; unknown
#'   fields remain unknown. No metadata or registration is inferred.
#' @param reference Optional SegmentationSet reference for existing validation.
#' @return SimpleList(geometry, coordinate_info, validation). Geometry is sf
#'   with exact cell_id; validation is the existing [validateGeometry()] report.
#' @details Source keys must be explicitly retained names, not automatic row
#'   numbers. No coordinate repair, CRS stripping, transformation, union,
#'   clipping or geometry dropping occurs. Census cells lacking geometry remain
#'   unresolved. Structural, identity and declaration errors stop extraction;
#'   per-cell eligibility failures are retained in the validation report.
#'   Only geometry interoperability is supported; this helper constructs no
#'   molecule universe. Both SpatialFeatureExperiment and sf are optional.
#' @examples
#' if (requireNamespace("SpatialFeatureExperiment", quietly = TRUE) &&
#'     requireNamespace("sf", quietly = TRUE)) {
#'     ring <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
#'     g <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(ring))))
#'     rownames(g) <- "source-cell"
#'     sfe <- SpatialFeatureExperiment::SpatialFeatureExperiment(
#'         assays = list(counts = matrix(0L, 1L, 1L,
#'             dimnames = list("g", "source-cell"))),
#'         colData = S4Vectors::DataFrame(sample_id = "roi"),
#'         spatialCoords = matrix(c(0.5, 0.5), 1L,
#'             dimnames = list(NULL, c("x", "y"))),
#'         sample_id = "roi", unit = "micron")
#'     SpatialFeatureExperiment::colGeometry(sfe, "cellSeg", sample_id = "roi",
#'         translate = FALSE) <- g
#'     sidecar <- geometryFromSpatialFeatureExperiment(sfe, "cellSeg", "roi",
#'         "col", data.frame(geometry_id = "source-cell", cell_id = "cell1"),
#'         data.frame(cell_id = c("cell1", "empty")),
#'         S4Vectors::SimpleList(units = "micron", specimen_id = "roi1"))
#'     sidecar$validation
#' }
#' @export
geometryFromSpatialFeatureExperiment <- function(object, type, sample_id,
    geometry_kind, cell_map, cells, coordinate_info, reference = NULL) {
    .adapterPackage("SpatialFeatureExperiment", object)
    if (!requireNamespace("sf", quietly = TRUE)) {
        stop("Optional package sf is required for geometry extraction",
             call. = FALSE)
    }
    geometry <- .sfeSourceGeometry(object, type, sample_id, geometry_kind)
    geometry <- .sfeMappedGeometry(geometry, cell_map)
    if (!methods::is(cells, "SegmentationPartition")) {
        cells <- .ingestionTable(cells, "cells")
    }
    report <- validateGeometry(geometry, cells, coordinate_info, reference)
    .adapterValidation(report)
    S4Vectors::SimpleList(geometry = geometry,
        coordinate_info = coordinate_info, validation = report)
}
