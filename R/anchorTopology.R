# Assisted-by: OpenAI Codex
.topologyRequire <- function(report) {
    failed <- report$status != "PASS"
    if (any(failed)) {
        stop("Topology-analysis prerequisite did not PASS: ",
            paste(paste0(report$check[failed], " [", report$status[failed],
                         "]: ", report$message[failed]), collapse = "; "),
            call. = FALSE)
    }
    report
}

.topologyWindow <- function(window) {
    if (!.baseMeasurement(window) || !is.null(dim(window)) ||
        length(window) != 4L || is.null(names(window)) ||
        anyNA(names(window)) || anyDuplicated(names(window)) ||
        !setequal(names(window), c("xmin", "xmax", "ymin", "ymax")) ||
        any(!is.finite(window)) || window["xmin"] >= window["xmax"] ||
        window["ymin"] >= window["ymax"]) {
        stop("window requires four finite named numeric bounds with ",
             "xmin < xmax and ymin < ymax", call. = FALSE)
    }
    window
}

.topologyCoordinates <- function(comparison, info) {
    stored <- comparison@coordinate_info
    if (is.null(info)) info <- stored
    declarations <- .coordinateDeclarations(info, stored = stored,
        specimen = comparison@metadata$specimen_id, prefix = "topology")
    .topologyRequire(declarations[declarations$status == "FAIL" &
                                      declarations$severity == "error", ])
    required <- .validationCheck("topology.coordinate_prerequisites", {
        shape <- .coordinateListProblems(info)
        if (!isTRUE(shape)) stop(shape)
        for (field in c("units", "origin", "axis_directions", "transforms",
                        "registration", "specimen_id")) {
            value <- info[[field]]
            known_none <- field == "transforms" &&
                .structuredValue(value) && !length(value)
            if (is.null(value) || (!known_none && .unknownMetadata(value))) {
                stop("Required common coordinate declaration unknown: ", field)
            }
        }
        TRUE
    })
    .topologyRequire(required)
    checked <- declarations[declarations$check !=
                                "topology.registration_evidence", ]
    # CRS is optional for explicitly declared local planar coordinates.
    checked <- checked[checked$check != "topology.crs", ]
    .topologyRequire(checked)
    list(info = info, report = .validationBind(list(declarations, required)))
}

.topologyCRS <- function(geometry, info) {
    .validationCheck("topology.crs.compatibility", {
        supplied <- geometry[!vapply(geometry, is.null, logical(1))]
        crs <- lapply(supplied, sf::st_crs)
        known <- vapply(crs, function(x) !is.na(x), logical(1))
        if (any(known)) {
            if (!all(known) || is.null(info$crs) ||
                .unknownMetadata(info$crs)) {
                stop("Known geometry CRS requires a common declared CRS")
            }
            if (!all(vapply(crs, function(x) isTRUE(x == crs[[1L]]),
                            logical(1)))) {
                stop("Supplied geometry CRS declarations conflict")
            }
        }
        if (any(vapply(supplied, function(g) isTRUE(sf::st_is_longlat(g)),
                       logical(1)))) {
            stop("Geographic geometry cannot PASS the planar area prerequisite")
        }
        TRUE
    })
}

.topologyGeometry <- function(geometry, cells, info, window) {
    report <- validateGeometry(geometry, cells, coordinate_info = info,
                               strict = TRUE)
    n <- nrow(cells)
    status <- rep("unavailable", n)
    reason <- rep("Census cell geometry unavailable", n)
    eligible <- interior <- rep(NA, n)
    present <- rep(FALSE, n)
    if (!is.null(geometry)) {
        index <- .matchIds(cells$cell_id, geometry$cell_id)
        present <- !is.na(index)
        for (i in which(present)) {
            prefix <- paste0("geometry.cell[", i, "].")
            checks <- report[startsWith(report$check, prefix), ]
            if (any(checks$status == "FAIL")) {
                status[i] <- "ineligible"
                eligible[i] <- FALSE
                reason[i] <- paste(checks$message[checks$status != "PASS"],
                                   collapse = "; ")
            } else if (!nrow(checks) || any(checks$status != "PASS")) {
                status[i] <- "unchecked"
                reason[i] <- paste(checks$message[checks$status != "PASS"],
                                   collapse = "; ")
            } else {
                bbox <- NULL
                box_check <- .validationCheck(paste0(prefix, "bbox"), {
                    bbox <- sf::st_bbox(sf::st_geometry(geometry)[index[i]])
                    if (any(!is.finite(bbox))) stop("Nonfinite polygon bbox")
                    TRUE
                }, severity = "warning")
                report <- .validationBind(list(report, box_check))
                if (.validationOK(box_check)) {
                    status[i] <- "eligible"
                    reason[i] <- "Geometry eligibility checks passed"
                    eligible[i] <- TRUE
                    interior[i] <- bbox["xmin"] > window["xmin"] &&
                        bbox["xmax"] < window["xmax"] &&
                        bbox["ymin"] > window["ymin"] &&
                        bbox["ymax"] < window["ymax"]
                } else {
                    status[i] <- "unchecked"
                    reason[i] <- box_check$message
                }
            }
        }
    }
    cells$geometry_present <- present
    cells$geometry_status <- status
    cells$geometry_reason <- reason
    cells$geometry_eligible <- eligible
    cells$strict_interior <- interior
    list(cells = cells, report = report)
}

.topologyStrata <- function(cells, index, other, other_index, reciprocal,
                            anchors, topology) {
    strata <- rep("non_reciprocal_matched", nrow(cells))
    strata[cells$dominance_count > 1L] <- "tied_dominance"
    strata[cells$partner_count == 0L] <- "zero_correspondence"
    for (e in which(reciprocal)) {
        i <- index[e]
        j <- other_index[e]
        statuses <- c(cells$geometry_status[i], other$geometry_status[j])
        strata[i] <- if ("unavailable" %in% statuses) {
            "reciprocal_geometry_unavailable"
        } else if ("ineligible" %in% statuses) {
            "reciprocal_geometry_ineligible"
        } else if ("unchecked" %in% statuses) {
            "reciprocal_geometry_unchecked"
        } else if (anchors[e]) {
            "interior_anchor"
        } else {
            "reciprocal_noninterior"
        }
    }
    cells$selection_stratum <- strata
    cells$interior_anchor <- rep(FALSE, nrow(cells))
    cells$interior_anchor[index[anchors]] <- TRUE
    cells$anchor_topology <- rep("undefined", nrow(cells))
    cells$anchor_topology[index[anchors]] <- topology[anchors]
    cells
}

.topologyResult <- function(comparison, cells, validation, window) {
    h <- comparison@edges
    q <- S4Vectors::from(h)
    s <- S4Vectors::to(h)
    c1 <- cells[[1L]]
    c2 <- cells[[2L]]
    eligible <- c1$geometry_eligible[q] & c2$geometry_eligible[s]
    interior <- eligible & c1$strict_interior[q] & c2$strict_interior[s]
    included <- !is.na(interior) & interior
    e <- S4Vectors::mcols(h)
    anchors <- e$reciprocal_dominant & included
    # Count unique reverse/forward dominance, including nonreciprocal relations.
    fan_out <- tabulate(q[included & e$unique_dominant_2], nrow(c1))
    fan_in <- tabulate(s[included & e$unique_dominant_1], nrow(c2))
    topology <- rep("undefined", length(h))
    topology[anchors] <- ifelse(fan_out[q[anchors]] > 1L,
        ifelse(fan_in[s[anchors]] > 1L, "complex", "split_like"),
        ifelse(fan_in[s[anchors]] > 1L, "merge_like", "one_to_one_like"))
    c1$fan_out <- fan_out
    c2$fan_in <- fan_in
    cells[[1L]] <- .topologyStrata(c1, q, c2, s, e$reciprocal_dominant,
                                  anchors, topology)
    cells[[2L]] <- .topologyStrata(c2, s, c1, q, e$reciprocal_dominant,
                                  anchors, topology)
    e$geometry_eligible <- eligible
    e$interior_edge <- interior
    e$interior_anchor <- anchors
    e$fan_out <- fan_out[q]
    e$fan_in <- fan_in[s]
    e$anchor_topology <- topology
    S4Vectors::mcols(h) <- e
    a <- S4Vectors::DataFrame(edge_index = which(anchors))
    a$cell_id_1 <- c1$cell_id[q[anchors]]
    a$cell_id_2 <- c2$cell_id[s[anchors]]
    a$fan_out <- fan_out[q[anchors]]
    a$fan_in <- fan_in[s[anchors]]
    a$anchor_topology <- topology[anchors]
    S4Vectors::SimpleList(anchors = a, cells = cells, edges = h,
                         validation = validation, window = window)
}

#' Describe secondary reciprocal-anchor topology
#'
#' Select unique reciprocal dominant edges of a complete molecule correspondence
#' graph whose endpoints have eligible geometry strictly inside a rectangular
#' observation window. Classes describe dominance relationships, not biological
#' truth, segmentation accuracy or method superiority.
#' @param comparison A SegmentationComparison with complete cell censuses.
#' @param geometry A named list with exactly both segmentation names. Each entry
#'   is an sf object with explicit typed cell_id, or NULL for unavailable
#'   geometry.
#' @param window Required finite unclassed numeric rectangle with exactly one
#'   each of xmin, xmax, ymin, ymax; xmin < xmax and ymin < ymax.
#' @param coordinate_info Optional common coordinate declaration SimpleList.
#'   NULL uses comparison declarations. Known compatible units, origin,
#'   axis_directions, transforms, registration and specimen_id are required.
#'   An empty structured transforms value explicitly declares no transforms.
#'   Local planar coordinates may have undeclared CRS; known geometry CRS must
#'   match a common declaration. Registration/source evidence remains unchecked.
#' @return A SimpleList with anchors (DataFrame), cells (named SimpleList of
#'   complete cell-summary DataFrames), edges (complete Hits with augmented
#'   metadata), validation (comparison/coordinates/named geometry reports) and
#'   window (the validated input). Anchors contain edge_index, typed cell_id_1
#'   and cell_id_2, integer fan_out and fan_in, and character anchor_topology.
#'   Cells add geometry_present, geometry_status, geometry_reason,
#'   geometry_eligible, strict_interior, selection_stratum, interior_anchor,
#'   anchor_topology and side-1 fan_out or side-2 fan_in. Edges add
#'   geometry_eligible, interior_edge, interior_anchor, fan_out, fan_in and
#'   anchor_topology. Original census and canonical edge order are preserved.
#' @details
#' Geometry validation is reused without repair. Missing geometry is
#'   unavailable;
#'   invalid, empty, nonfinite or nonpositive-area geometry is ineligible;
#'   UNKNOWN backend prerequisites remain unchecked. These cells and all edges
#'   are retained. Geometry properties not checked remain NA, and only proven
#'   eligible interior cells contribute to anchors or fans. Unsupported types,
#'   malformed or unknown IDs and known declaration conflicts stop analysis.
#'
#' Strict interior compares the complete polygon bbox to all four bounds using
#'   strict inequalities. Boundary-touching, crossing and outside polygons are
#'   non-interior. Holes and all multipolygon parts are retained. Interior
#'   status
#'   does not certify independence from unavailable exterior molecule context.
#'
#' Fan-out counts eligible interior side-2 cells whose unique complete-graph
#'   dominant side-1 partner is the anchor's side-1 cell. Fan-in reverses these
#'   roles. Both endpoints of each contributing relation must be interior, but
#'   contributing relations need not be reciprocal. Classes are exactly
#'   one_to_one_like (1,1), split_like (>1,1), merge_like (1,>1),
#'   complex (>1,>1).
#'   Non-anchor topology is undefined. Selection strata distinguish isolates,
#'   ties, nonreciprocity, reciprocal unavailable/ineligible/unchecked geometry,
#'   reciprocal non-interior pairs and interior anchors. No adjacency or
#'   distance
#'   predicate, support filter or state input changes the complete primary
#'   graph.
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'     m <- S4Vectors::DataFrame(molecule_id = "m", feature_id = "g")
#'     p <- SegmentationPartition(m, "first",
#'         S4Vectors::DataFrame(cell_id = "x"), "m", "assigned", "x")
#'     q <- SegmentationPartition(m, "second",
#'         S4Vectors::DataFrame(cell_id = "a"), "m", "assigned", "a")
#'     ss <- SegmentationSet(m, list(first = p, second = q), "specimen")
#'     cmp <- compareSegmentations(ss, "first", "second")
#'     ring <- rbind(c(1, 1), c(2, 1), c(2, 2), c(1, 2), c(1, 1))
#'     shape <- sf::st_sfc(sf::st_polygon(list(ring)))
#'     geometry <- list(first = sf::st_sf(cell_id = "x", geometry = shape),
#'                      second = sf::st_sf(cell_id = "a", geometry = shape))
#'     info <- S4Vectors::SimpleList(units = "micrometres", origin = c(0, 0),
#'         axis_directions = c("right", "up"),
#'         transforms = S4Vectors::SimpleList(), registration = "common XY",
#'         specimen_id = "specimen")
#'     anchorTopology(cmp, geometry,
#'         c(xmin = 0, xmax = 10, ymin = 0, ymax = 10), info)$anchors
#' }
#' @export
anchorTopology <- function(comparison, geometry, window,
                           coordinate_info = NULL) {
    validation <- .stateComparisonChecks(comparison, prefix = "topology",
                                         enforce = .topologyRequire)
    window <- .topologyWindow(window)
    if (!.geometryAvailable()) {
        stop("anchorTopology requires optional sf", call. = FALSE)
    }
    sides <- comparison@segmentation_names
    if (!is.list(geometry) || is.object(geometry) ||
        is.null(names(geometry)) || anyNA(names(geometry)) ||
        anyDuplicated(names(geometry)) || length(geometry) != 2L ||
        !setequal(names(geometry), sides)) {
        stop("geometry must be a named list of both segmentations",
             call. = FALSE)
    }
    geometry <- geometry[sides]
    coordinates <- .topologyCoordinates(comparison, coordinate_info)
    checked <- lapply(seq_along(sides), function(i) {
        cells <- if (i == 1L) comparison@cell_summary_1 else
            comparison@cell_summary_2
        .topologyGeometry(geometry[[i]], cells, coordinates$info, window)
    })
    crs <- .topologyCRS(geometry, coordinates$info)
    .topologyRequire(crs)
    cells <- S4Vectors::SimpleList(lapply(checked, `[[`, "cells"))
    reports <- S4Vectors::SimpleList(lapply(checked, `[[`, "report"))
    names(cells) <- names(reports) <- sides
    validation <- S4Vectors::SimpleList(comparison = .validationBind(
        list(validation, .validationUnknown("topology.source_evidence",
            "Source identity/census/authenticity not reverified"))),
        coordinates = .validationBind(list(coordinates$report, crs)),
        geometry = reports)
    .topologyResult(comparison, cells, validation, window)
}
