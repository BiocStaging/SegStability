.geometryRequire <- function(report) {
    failed <- report$status != "PASS"
    if (any(failed)) {
        stop("Geometry-comparison prerequisite did not PASS: ",
            paste(paste0(report$check[failed], " [", report$status[failed],
                         "]: ", report$message[failed]), collapse = "; "),
            call. = FALSE)
    }
    report
}

.geometryCensus <- function(geometry, cells, info) {
    report <- validateGeometry(geometry, cells, coordinate_info = info,
                               strict = TRUE)
    n <- nrow(cells)
    index <- rep(NA_integer_, n)
    if (!is.null(geometry)) index <- .matchIds(cells$cell_id, geometry$cell_id)
    status <- rep("unavailable", n)
    reason <- rep("Census cell geometry unavailable", n)
    eligible <- rep(NA, n)
    area <- rep(NA_real_, n)
    for (i in which(!is.na(index))) {
        prefix <- paste0("geometry.cell[", i, "].")
        checks <- report[startsWith(report$check, prefix), ]
        if (any(checks$status == "FAIL")) {
            status[i] <- "ineligible"
            eligible[i] <- FALSE
        } else if (!nrow(checks) || any(checks$status != "PASS")) {
            status[i] <- "unchecked"
        } else {
            value <- NA_real_
            measured <- .validationCheck(paste0(prefix, "measurement"), {
                value <- as.numeric(sf::st_area(
                    sf::st_geometry(geometry)[index[i]]))
                is.finite(value) && value > 0
            }, severity = "warning")
            report <- .validationBind(list(report, measured))
            checks <- rbind(checks, measured)
            if (.validationOK(measured)) {
                status[i] <- "eligible"
                eligible[i] <- TRUE
                area[i] <- value
            } else {
                status[i] <- "unchecked"
            }
        }
        reason[i] <- if (status[i] == "eligible") {
            "Geometry eligibility checks passed"
        } else {
            paste(checks$message[checks$status != "PASS"], collapse = "; ")
        }
    }
    cells$geometry_present <- !is.na(index)
    cells$geometry_status <- status
    cells$geometry_reason <- reason
    cells$geometry_eligible <- eligible
    cells$area <- area
    list(cells = cells, index = index, report = report)
}

.geometryCandidates <- function(geometry, checked) {
    q <- which(checked[[1L]]$cells$geometry_status == "eligible")
    s <- which(checked[[2L]]$cells$geometry_status == "eligible")
    pairs <- list(q = integer(), s = integer())
    report <- .validationCheck("geometry_comparison.enumeration", {
        if (length(q) && length(s)) {
            left <- sf::st_geometry(geometry[[1L]])[checked[[1L]]$index[q]]
            right <- sf::st_geometry(geometry[[2L]])[checked[[2L]]$index[s]]
            hits <- sf::st_intersects(left, right, sparse = TRUE)
            pairs$q <- rep(q, lengths(hits))
            pairs$s <- s[unlist(hits, use.names = FALSE)]
        }
        TRUE
    })
    # Unknown enumeration cannot certify the all-positive-pair population.
    .geometryRequire(report)
    list(pairs = pairs, report = report)
}

.geometryPairTable <- function(comparison, candidates, cells) {
    h <- comparison@edges
    q <- c(S4Vectors::from(h), candidates$q)
    s <- c(S4Vectors::to(h), candidates$s)
    order <- order(q, s, method = "radix")
    q <- q[order]
    s <- s[order]
    key <- paste(q, s, sep = ":")
    keep <- !duplicated(key)
    q <- q[keep]
    s <- s[keep]
    edge_key <- paste(S4Vectors::from(h), S4Vectors::to(h), sep = ":")
    out <- S4Vectors::DataFrame(cell_index_1 = q, cell_index_2 = s,
        cell_id_1 = cells[[1L]]$cell_id[q], cell_id_2 = cells[[2L]]$cell_id[s],
        component_id_1 = cells[[1L]]$component_id[q],
        component_id_2 = cells[[2L]]$component_id[s],
        edge_index = match(key[keep], edge_key))
    out$molecule_edge <- !is.na(out$edge_index)
    out
}

.geometryMeasurePair <- function(left, right, a, b) {
    intersection <- sf::st_intersection(left, right)
    shared <- sum(as.numeric(sf::st_area(intersection)))
    union <- max(a, b) + (min(a, b) - shared)
    if (!is.finite(shared) || shared < 0 || shared > min(a, b) ||
        !is.finite(union) || union <= 0) {
        stop("Intersection/union area cannot support safe overlap metrics")
    }
    # Scale Dice to avoid overflowing 2 * intersection or area_1 + area_2.
    scale <- max(a, b)
    c(a, b, shared, union, shared / union,
      2 * (shared / scale) / (a / scale + b / scale), shared / a, shared / b)
}

.geometryMetrics <- function(pairs, geometry, checked) {
    n <- nrow(pairs)
    fields <- c("area_1", "area_2", "intersection_area", "union_area",
                "iou", "dice", "overlap_fraction_1", "overlap_fraction_2")
    for (field in fields) pairs[[field]] <- rep(NA_real_, n)
    status <- reason <- rep("", n)
    reports <- list()
    for (e in seq_len(n)) {
        q <- pairs$cell_index_1[e]
        s <- pairs$cell_index_2[e]
        statuses <- c(checked[[1L]]$cells$geometry_status[q],
                      checked[[2L]]$cells$geometry_status[s])
        if (!all(statuses == "eligible")) {
            status[e] <- if ("unavailable" %in% statuses) "unavailable" else
                if ("ineligible" %in% statuses) "ineligible" else "unchecked"
            reason[e] <- paste(checked[[1L]]$cells$geometry_reason[q],
                               checked[[2L]]$cells$geometry_reason[s],
                               sep = "; ")
            next
        }
        values <- NULL
        label <- paste0("geometry_comparison.pair[", e, "]")
        report <- .validationCheck(label, {
            values <- .geometryMeasurePair(
                sf::st_geometry(geometry[[1L]])[checked[[1L]]$index[q]],
                sf::st_geometry(geometry[[2L]])[checked[[2L]]$index[s]],
                checked[[1L]]$cells$area[q], checked[[2L]]$cells$area[s])
            TRUE
        }, severity = "warning")
        reports[[length(reports) + 1L]] <- report
        status[e] <- if (.validationOK(report)) "eligible" else "unchecked"
        reason[e] <- report$message
        if (.validationOK(report)) {
            for (i in seq_along(fields)) pairs[[fields[i]]][e] <- values[i]
        }
    }
    pairs$geometry_status <- status
    pairs$geometry_reason <- reason
    pairs$positive_overlap <- ifelse(status == "eligible",
                                    pairs$intersection_area > 0, NA)
    # Known zero-area nonedges are outside the positive-overlap population.
    # Failed candidate measurements stay explicit.
    keep <- pairs$molecule_edge | is.na(pairs$positive_overlap) |
        pairs$positive_overlap
    pairs <- pairs[keep, ]
    pairs$pair_index <- seq_len(nrow(pairs))
    report <- if (length(reports)) .validationBind(reports) else
        .validationReport("geometry_comparison.metrics", "NOT_APPLICABLE",
                          message = "No eligible pairs to measure")
    list(pairs = pairs, report = report)
}

.geometryResult <- function(comparison, pairs, cells, validation, info) {
    computable <- !is.na(pairs$positive_overlap) & pairs$positive_overlap
    complete_metrics <- !any(pairs$geometry_status == "unchecked") &&
        !any(vapply(cells, function(x) any(x$geometry_status == "unchecked"),
                    logical(1)))
    selected <- integer()
    if (complete_metrics) {
        eligible <- which(computable)
        matched <- .geometryMatching(pairs$cell_index_1[eligible],
            pairs$cell_index_2[eligible], pairs$iou[eligible])
        selected <- eligible[matched]
    }
    pairs$geometry_matched <- rep(if (complete_metrics) FALSE else NA,
                                 nrow(pairs))
    pairs$geometry_matched[selected] <- TRUE
    for (side in 1:2) {
        matched <- rep(NA, nrow(cells[[side]]))
        eligible <- cells[[side]]$geometry_status == "eligible"
        if (complete_metrics) matched[eligible] <- FALSE
        matched[pairs[[paste0("cell_index_", side)]][selected]] <- TRUE
        cells[[side]]$geometry_matched <- matched
    }
    h <- comparison@edges
    index <- match(seq_along(h), pairs$edge_index)
    e <- S4Vectors::mcols(h)
    fields <- c("pair_index", "geometry_status", "geometry_reason",
                "positive_overlap", "geometry_matched", "area_1", "area_2",
                "intersection_area", "union_area", "iou", "dice",
                "overlap_fraction_1", "overlap_fraction_2")
    for (field in fields) e[[field]] <- pairs[[field]][index]
    S4Vectors::mcols(h) <- e
    coverage <- all(vapply(cells, function(x)
        all(x$geometry_status == "eligible"), logical(1)))
    summary <- S4Vectors::SimpleList(coverage_complete = coverage,
        matching_status = if (!complete_metrics) "unchecked" else
            if (coverage) "complete" else "eligible_subset",
        n_molecule_edges = length(h),
        n_positive_geometry_pairs = sum(computable),
        n_geometry_only_pairs = sum(computable & !pairs$molecule_edge),
        n_matches = if (complete_metrics) length(selected) else NA_integer_,
        matching_objective = if (complete_metrics) sum(pairs$iou[selected]) else
            NA_real_,
        n_geometry_only_matches = if (complete_metrics)
            sum(!pairs$molecule_edge[selected]) else NA_integer_,
        n_molecule_edges_omitted = if (complete_metrics)
            length(h) - sum(pairs$molecule_edge[selected]) else NA_integer_)
    statuses <- c("eligible", "unavailable", "ineligible", "unchecked")
    counts <- lapply(cells, function(x) S4Vectors::DataFrame(
        geometry_status = statuses,
        n_cells = unname(vapply(statuses,
            function(status) sum(x$geometry_status == status), integer(1)))))
    summary$cells <- S4Vectors::SimpleList(counts)
    validation$matching <- if (complete_metrics) {
        .validationReport("geometry_comparison.matching", "PASS",
            message = "Exact maximum-IoU objective on eligible positive pairs")
    } else {
        .validationUnknown("geometry_comparison.matching",
            "Unknown geometry/metric checks preclude certified matching")
    }
    S4Vectors::SimpleList(pairs = pairs, matching = pairs[selected, ],
        cells = cells, edges = h, components = comparison@components,
        summary = summary, validation = validation, coordinate_info = info)
}

#' Describe polygon agreement and derived maximum-IoU matching
#'
#' Compare supplied polygons across the full molecule correspondence and all
#' positive-area geometric pairs. Geometry matching is a descriptive annotation;
#' it does not replace molecule correspondence or rank segmentation accuracy.
#' @param comparison A SegmentationComparison with complete cell censuses.
#' @param geometry A base named list with exactly both segmentation names, each
#'   containing an sf object with explicit typed cell_id, or NULL.
#' @param coordinate_info Common coordinate declaration SimpleList. NULL uses
#'   comparison declarations. Known compatible units, origin, axis_directions,
#'   transforms, registration and specimen_id are required. Local planar CRS may
#'   be undeclared; known CRS must match a common declaration. Geographic
#'   coordinates are ineligible for these planar measurements.
#' @return A SimpleList with pairs and matching DataFrames, named full-census
#'   cells SimpleList, original edges Hits with derived annotations, unchanged
#'   components, summary SimpleList, validation reports and coordinate_info.
#'   Pairs contain pair_index, cell_index_1/2, typed cell_id_1/2, original
#'   component_id_1/2, nullable edge_index, molecule_edge, positive_overlap,
#'   geometry_matched, geometry_status/reason and the eight metric fields below.
#'   Matching contains the selected pair rows. Cells add geometry_present,
#'   geometry_status/reason, geometry_eligible, area and geometry_matched.
#'   Edges add pair_index, geometry status/reason, overlap/matching flags and
#'   metrics, preserving every original field and graph relation.
#' @details
#' area_1 and area_2 are polygon areas; intersection_area is their intersection;
#' union_area is their area union. Areas use declared units squared. iou (also
#' Jaccard) is intersection/union; dice is twice intersection/summed areas;
#' overlap_fraction_1/2 divide intersection by the respective endpoint area.
#' Ratios are dimensionless. Eligible disjoint molecule edges have zero overlap.
#' Nonpositive areas never enter metrics. Missing, invalid, empty, nonfinite and
#' unchecked geometry remains explicit; unresolved pair metrics are NA.
#'
#' Pairs retain every molecule edge and all eligible positive-area geometry
#' pairs, including those without shared molecules. Failed measurements remain
#' explicit. No window or interior filter applies. Split, merge and complex
#' correspondence stays pairwise; holes and all multipolygon parts are retained.
#'
#' Matching maximizes total unmodified IoU on eligible positive geometry pairs.
#' Error-free sums distinguish exact ties from near ties in supplied doubles.
#' Equal optimal objectives select the lexicographically first sorted list of
#' census-index pairs. No epsilon, weight perturbation or zero-weight match is
#' used. Missing/ineligible geometry gives matching_status eligible_subset;
#' unchecked operations give unchecked with no certified matching. Empty known
#' graphs have objective zero. Source registration/authenticity stays UNKNOWN.
#'
#' Summary reports coverage_complete, matching_status, n_molecule_edges,
#' n_positive_geometry_pairs, n_geometry_only_pairs, n_matches,
#' matching_objective, n_geometry_only_matches, n_molecule_edges_omitted and
#' named per-side cells status counts. The displayed objective is rounded to a
#' double; optimization and equality decisions use error-free sums internally.
#' Enumeration failure stops rather than claiming all positive pairs were found.
#' No repair, buffering, clipping, simplification, transformation or inferred
#' geometry/registration is performed. Optional sf is required only here.
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'     m <- S4Vectors::DataFrame(molecule_id = "m", feature_id = "g")
#'     p <- SegmentationPartition(m, "first",
#'         S4Vectors::DataFrame(cell_id = "x"), "m", "assigned", "x")
#'     q <- SegmentationPartition(m, "second",
#'         S4Vectors::DataFrame(cell_id = "a"), "m", "assigned", "a")
#'     cmp <- compareSegmentations(
#'         SegmentationSet(m, list(first = p, second = q), "specimen"),
#'         "first", "second")
#'     ring <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
#'     shape <- sf::st_sfc(sf::st_polygon(list(ring)))
#'     g <- list(first = sf::st_sf(cell_id = "x", geometry = shape),
#'               second = sf::st_sf(cell_id = "a", geometry = shape))
#'     info <- S4Vectors::SimpleList(units = "micrometres", origin = c(0, 0),
#'         axis_directions = c("right", "up"),
#'         transforms = S4Vectors::SimpleList(), registration = "common XY",
#'         specimen_id = "specimen")
#'     compareGeometry(cmp, g, info)$pairs
#' }
#' @export
compareGeometry <- function(comparison, geometry, coordinate_info = NULL) {
    checks <- .stateComparisonChecks(comparison, prefix = "geometry_comparison",
                                     enforce = .geometryRequire)
    if (!.geometryAvailable()) {
        stop("compareGeometry requires optional sf", call. = FALSE)
    }
    sides <- comparison@segmentation_names
    if (!is.list(geometry) || is.object(geometry) || length(geometry) != 2L ||
        is.null(names(geometry)) || anyNA(names(geometry)) ||
        anyDuplicated(names(geometry)) || !setequal(names(geometry), sides)) {
        stop("geometry must be a named list of both segmentations",
             call. = FALSE)
    }
    geometry <- geometry[sides]
    coordinates <- .topologyCoordinates(comparison, coordinate_info)
    checked <- lapply(1:2, function(i) .geometryCensus(geometry[[i]],
        if (i == 1L) comparison@cell_summary_1 else comparison@cell_summary_2,
        coordinates$info))
    crs <- .topologyCRS(geometry, coordinates$info)
    .geometryRequire(crs)
    cells <- S4Vectors::SimpleList(lapply(checked, `[[`, "cells"))
    reports <- S4Vectors::SimpleList(lapply(checked, `[[`, "report"))
    names(cells) <- names(reports) <- sides
    candidates <- .geometryCandidates(geometry, checked)
    pairs <- .geometryPairTable(comparison, candidates$pairs, cells)
    measured <- .geometryMetrics(pairs, geometry, checked)
    validation <- S4Vectors::SimpleList(comparison = .validationBind(
        list(checks, .validationUnknown("geometry_comparison.source_evidence",
            "Source identity/census/authenticity not reverified"))),
        coordinates = .validationBind(list(coordinates$report, crs)),
        geometry = reports, enumeration = candidates$report,
        metrics = measured$report)
    .geometryResult(comparison, measured$pairs, cells, validation,
                    coordinates$info)
}
