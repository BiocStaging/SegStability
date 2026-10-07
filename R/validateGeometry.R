.geometryAvailable <- function() requireNamespace("sf", quietly = TRUE)

.geometryCellChecks <- function(geometry, idx, check) {
    sfc <- sf::st_geometry(geometry)[idx]
    finite <- .validationCheck(paste0(check, ".finite"), {
        coordinates <- sf::st_coordinates(sfc)
        if (any(!is.finite(coordinates))) {
            stop("Nonfinite geometry coordinates; geometry is ineligible")
        }
        TRUE
    }, severity = "warning")
    empty <- .validationCheck(paste0(check, ".nonempty"), {
        if (isTRUE(sf::st_is_empty(sfc))) {
            stop("Empty geometry is ineligible")
        }
        TRUE
    }, severity = "warning")
    out <- list(finite, empty)
    if (.validationOK(finite) && .validationOK(empty)) {
        valid <- .validationCheck(paste0(check, ".valid"), {
            reason <- sf::st_is_valid(sfc, reason = TRUE)
            if (length(reason) != 1L || is.na(reason) ||
                !identical(reason, "Valid Geometry")) {
                stop("Invalid geometry is ineligible: ", reason)
            }
            TRUE
        }, severity = "warning")
        out[[length(out) + 1L]] <- valid
        if (.validationOK(valid)) {
            if (isTRUE(sf::st_is_longlat(sfc))) {
                out[[length(out) + 1L]] <- .validationUnknown(
                    paste0(check, ".area"),
                    "Planar area gate unavailable for geographic coordinates")
            } else {
                out[[length(out) + 1L]] <- .validationCheck(
                    paste0(check, ".area"), {
                        area <- as.numeric(sf::st_area(sfc))
                        if (length(area) != 1L || !is.finite(area) ||
                            area <= 0) {
                            stop("Nonpositive/nonfinite area is ineligible")
                        }
                        TRUE
                    }, severity = "warning")
            }
        } else {
            out[[length(out) + 1L]] <- .validationUnknown(
                paste0(check, ".area"),
                "Validity gate did not pass; area unchecked")
        }
    } else {
        out[[length(out) + 1L]] <- .validationUnknown(
            paste0(check, ".valid"),
            "Validity/area unavailable for empty or nonfinite geometry")
    }
    .validationBind(out)
}

.geometryCRSChecks <- function(geometry, info) {
    out <- list()
    if (.namedValidationList(info) && !is.null(info$crs) &&
        .coordinateShape("crs", info$crs) && !.unknownMetadata(info$crs)) {
        out[[1L]] <- .validationCheck("geometry.crs.object", {
            declared <- sf::st_crs(info$crs)
            actual <- sf::st_crs(geometry)
            if (is.na(declared) || is.na(actual) ||
                !isTRUE(actual == declared)) {
                stop("Geometry CRS conflicts with its declared CRS")
            }
            TRUE
        })
    }
    if (!length(out)) {
        out[[1L]] <- .validationUnknown("geometry.crs.object",
            "CRS linkage not checked without a known declaration")
    }
    if (isTRUE(sf::st_is_longlat(geometry))) {
        out[[length(out) + 1L]] <- .validationCheck("geometry.crs.local", {
            if (.namedValidationList(info) && is.character(info$units) &&
                any(info$units %in% c("micron", "micrometre", "micrometres",
                    "micrometer", "micrometers", "um", "pixel", "pixels",
                    "mm", "metre", "metres", "meter", "meters", "m"))) {
                stop("Geographic CRS conflicts with declared local units")
            }
            TRUE
        })
    }
    .validationBind(out)
}

#' Validate optional cell geometry without repair
#'
#' Geometry is external to core objects. Polygons, multipolygons and holes are
#' supported; overlaps are allowed. Cell IDs are matched exactly to the complete
#' census. No geometry is repaired, buffered, clipped, simplified or
#'   transformed.
#' @param geometry An sf object with explicit cell_id, or NULL when absent.
#' @param cells A DataFrame containing cell_id or a SegmentationPartition
#'   census.
#' @param coordinate_info Optional uniquely named coordinate declaration
#'   SimpleList
#'   with the fields described in validateCoordinates().
#' @param reference Optional SegmentationSet for specimen/declaration
#'   consistency.
#' @param strict Logical scalar. Stop ONLY for FAIL/error findings. Invalid,
#'   empty, nonpositive-area and nonfinite supplied geometry is FAIL/warning
#'   and ineligible; these findings alone never cause strict mode to throw.
#' @return A DataFrame with character check, status, severity and message
#'   columns.
#'   Missing census geometry is UNKNOWN/note; absent geometry
#'   NOT_APPLICABLE/note.
#'   Supplied geometry without optional sf is unchecked UNKNOWN/warning.
#'   Unknown/duplicate IDs, unsupported types and known conflicts are
#'   FAIL/error.
#'   Geographic planar area is unchecked. Declaration consistency alone cannot
#'   attest source linkage or registration.
#' @details
#' Per-cell check labels geometry.cell[i] refer to original census row i.
#' Geometry rows are aligned by exact IDs, including cells without geometry.
#' These diagnostics do not return altered geometry or certify spatial sources.
#' @examples
#' validateGeometry(NULL, S4Vectors::DataFrame(cell_id = "c1"), strict = TRUE)
#' if (requireNamespace("sf", quietly = TRUE)) {
#'     ring <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
#'     g <- sf::st_sf(cell_id = "c1",
#'         geometry = sf::st_sfc(sf::st_polygon(list(ring))))
#'     validateGeometry(g, S4Vectors::DataFrame(cell_id = "c1"), strict = TRUE)
#' }
#' @export
validateGeometry <- function(geometry, cells, coordinate_info = NULL,
                             reference = NULL, strict = FALSE) {
    if (is.null(geometry)) {
        return(.validationFinish(.validationReport("geometry.presence",
            "NOT_APPLICABLE", message = "Optional geometry not supplied"),
            strict))
    }
    input <- .validationCheck("geometry.input",
        if (inherits(geometry, "sf")) TRUE else "Expected sf with cell_id")
    if (!.validationOK(input)) {
        return(.validationFinish(input, strict))
    }
    if (!.geometryAvailable()) {
        dependency <- .validationReport("geometry.dependency", "UNKNOWN",
            "warning", "Optional sf unavailable; geometry was not checked")
        return(.validationFinish(.validationBind(list(input, dependency)),
                                  strict))
    }
    structure <- .validationCheck("geometry.structure", {
        sfc <- sf::st_geometry(geometry)
        if (!inherits(sfc, "sfc") || length(sfc) != nrow(geometry)) {
            stop("Malformed sf geometry column or row count")
        }
        TRUE
    })
    if (!.validationOK(structure)) {
        return(.validationFinish(.validationBind(list(input, structure,
            .validationUnknown("geometry.scientific",
                "Geometry checks unavailable after structural failure"))),
            strict))
    }
    if (methods::is(cells, "SegmentationPartition")) {
        cells <- cells@cells
    }
    census <- .validationCheck("geometry.census", {
        if (!methods::is(cells, "DataFrame") ||
            !"cell_id" %in% colnames(cells) || anyDuplicated(colnames(cells))) {
            stop("Expected a census DataFrame with cell_id")
        }
        .validateIds(cells$cell_id, "cell census")
        TRUE
    })
    identity <- .validationCheck("geometry.identity", {
        if (!"cell_id" %in% colnames(geometry) ||
            anyDuplicated(colnames(geometry))) {
            stop("Geometry requires unique cell_id column")
        }
        .validateIds(geometry$cell_id, "geometry cell_id")
        TRUE
    })
    types <- .validationCheck("geometry.types", {
        if (any(!as.character(sf::st_geometry_type(geometry)) %in%
                c("POLYGON", "MULTIPOLYGON"))) {
            stop("Only polygons and multipolygons are supported")
        }
        TRUE
    })
    out <- list(input, structure, census, identity, types)
    eligible_checks <- FALSE
    if (.validationOK(census) && .validationOK(identity)) {
        membership <- .validationCheck("geometry.membership", {
            idx <- .matchIds(geometry$cell_id, cells$cell_id)
            if (anyNA(idx)) {
                stop("Unknown geometry cell IDs outside the census")
            }
            TRUE
        })
        out[[length(out) + 1L]] <- membership
        if (.validationOK(membership) && .validationOK(types)) {
            eligible_checks <- TRUE
            idx <- .matchIds(cells$cell_id, geometry$cell_id)
            for (i in seq_along(idx)) {
                check <- paste0("geometry.cell[", i, "]")
                out[[length(out) + 1L]] <- if (is.na(idx[i])) {
                    .validationUnknown(paste0(check, ".presence"),
                                       "Census cell geometry unavailable")
                } else {
                    .geometryCellChecks(geometry, idx[i], check)
                }
            }
        }
    }
    if (!eligible_checks) {
        out[[length(out) + 1L]] <- .validationUnknown("geometry.eligibility",
            "Cell eligibility unchecked after census/identity/type failure")
    }
    ref_info <- NULL
    if (!is.null(reference)) {
        ref <- .coordinateReference(reference)
        ref$check <- "geometry.reference"
        out[[length(out) + 1L]] <- ref
        if (.validationOK(ref)) {
            ref_info <- reference@coordinate_info
            ref_info$specimen_id <- reference@specimen_id
        }
    }
    out[[length(out) + 1L]] <- .coordinateDeclarations(coordinate_info,
        reference = ref_info, prefix = "geometry")
    out[[length(out) + 1L]] <- .geometryCRSChecks(geometry, coordinate_info)
    .validationFinish(.validationBind(out), strict)
}
