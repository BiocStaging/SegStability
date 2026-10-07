.exactCounts <- function(x) {
    (is.integer(x) || identical(class(x), "numeric")) &&
        !anyNA(x) && all(is.finite(x) & x >= 0 & x <= 2^53 & x == floor(x))
}

.countTotal <- function(x) {
    if (!.exactCounts(x)) return(NA_real_)
    total <- 0
    for (value in x) {
        # Guard before addition: 2^53 + 1 can round back to 2^53.
        if (value > 2^53 - total) return(NA_real_)
        total <- total + value
    }
    total
}

.sumByIndex <- function(index, mass, n) {
    result <- numeric(n)
    if (length(index)) {
        sums <- rowsum(mass, index, reorder = FALSE)
        result[as.integer(rownames(sums))] <- sums[, 1L]
    }
    result
}

.cellSummaryProblems <- function(cells) {
    required <- c("cell_id", "assigned_mass", "matched_mass", "unmatched_mass",
                  "partner_count", "dominance_count", "component_id",
                  "matched_fraction")
    if (!all(required %in% colnames(cells))) {
        return("Missing cell summary columns")
    }
    ids <- tryCatch({
        .validateIds(cells$cell_id, "cell summary")
        TRUE
    }, error = function(e) conditionMessage(e))
    if (!isTRUE(ids)) return(ids)
    counts <- required[2:6]
    if (!all(vapply(counts, function(nm) .exactCounts(cells[[nm]]),
                    logical(1))) ||
        any(cells$matched_mass > cells$assigned_mass) ||
        any(cells$unmatched_mass != cells$assigned_mass - cells$matched_mass) ||
        any(cells$dominance_count > cells$partner_count)) {
        return("Invalid cell masses or partner/dominance counts")
    }
    expected <- rep(NA_real_, nrow(cells))
    positive <- cells$assigned_mass > 0
    expected[positive] <- cells$matched_mass[positive] /
        cells$assigned_mass[positive]
    if (!is.numeric(cells$matched_fraction) ||
        !identical(as.numeric(cells$matched_fraction), expected)) {
        return("Invalid matched_fraction or zero-denominator behavior")
    }
    TRUE
}

.componentClass <- function(n1, n2) {
    # API section 8 names; scientific section 3 cardinality definitions.
    as.character(ifelse(n1 == 0L | n2 == 0L, "unmatched",
        ifelse(n1 == 1L & n2 == 1L, "one_to_one",
            ifelse(n1 == 1L, "split",
                ifelse(n2 == 1L, "merge", "many_to_many")))))
}

.componentProblems <- function(components, cells1, cells2) {
    required <- c("component_id", "n_segmentation_1", "n_segmentation_2",
                  "component_class")
    if (!all(required %in% colnames(components)) ||
        !identical(components$component_id, seq_len(nrow(components)))) {
        return("Invalid component table or component IDs")
    }
    for (cells in list(cells1, cells2)) {
        id <- cells$component_id
        if (!is.integer(id) || anyNA(id) ||
            any(id < 1L | id > nrow(components))) {
            return("Invalid cell component references")
        }
    }
    n1 <- tabulate(cells1$component_id, nrow(components))
    n2 <- tabulate(cells2$component_id, nrow(components))
    if (!.exactCounts(components$n_segmentation_1) ||
        !.exactCounts(components$n_segmentation_2) ||
        any(components$n_segmentation_1 != n1) ||
        any(components$n_segmentation_2 != n2) || any(n1 + n2 == 0L) ||
        any((n1 == 0L | n2 == 0L) & n1 + n2 != 1L) ||
        !identical(components$component_class, .componentClass(n1, n2))) {
        return("Invalid component cardinalities or component_class")
    }
    TRUE
}

.graphStructureProblems <- function(h, c1, c2) {
    q <- S4Vectors::from(h)
    s <- S4Vectors::to(h)
    if (S4Vectors::nLnode(h) != nrow(c1) ||
        S4Vectors::nRnode(h) != nrow(c2) || anyNA(q) || anyNA(s) ||
        any(q < 1L | q > nrow(c1)) || any(s < 1L | s > nrow(c2)) ||
        !identical(order(q, s, method = "radix"), seq_along(q)) ||
        any(q[-1L] == q[-length(q)] & s[-1L] == s[-length(s)])) {
        return("Invalid edge dimensions, references or canonical unique order")
    }
    e <- S4Vectors::mcols(h)
    flags <- c("dominant_1", "dominant_2", "unique_dominant_1",
               "unique_dominant_2", "reciprocal_dominant")
    if (!all(c("weight", "fraction_1", "fraction_2", flags) %in%
             colnames(e)) || is.na(.countTotal(e$weight)) ||
        any(e$weight == 0) ||
        !all(vapply(flags, function(n) is.logical(e[[n]]) &&
                    !anyNA(e[[n]]), logical(1)))) {
        return("Invalid positive whole edge weights or graph metadata")
    }
    TRUE
}

.directionalEdgeProblems <- function(h, c1, c2) {
    q <- S4Vectors::from(h)
    s <- S4Vectors::to(h)
    e <- S4Vectors::mcols(h)
    for (side in 1:2) {
        cells <- if (side == 1L) c1 else c2
        idx <- if (side == 1L) q else s
        f <- e[[paste0("fraction_", side)]]
        d <- e[[paste0("dominant_", side)]]
        u <- e[[paste0("unique_dominant_", side)]]
        if (any(cells$assigned_mass[idx] == 0) ||
            !is.numeric(f) || anyNA(f) || any(f < 0 | f > 1) ||
            any(f != e$weight / cells$assigned_mass[idx]) ||
            any(cells$matched_mass !=
                .sumByIndex(idx, e$weight, nrow(cells))) ||
            any(cells$partner_count != tabulate(idx, nrow(cells))) ||
            any(cells$dominance_count != tabulate(idx[d], nrow(cells))) ||
            any(u != (d & cells$dominance_count[idx] == 1L))) {
            return("Invalid edge fractions or cell correspondence accounting")
        }
    }
    if (any(e$reciprocal_dominant !=
            (e$unique_dominant_1 & e$unique_dominant_2))) {
        return("Invalid reciprocal dominance flags")
    }
    TRUE
}

.assignmentSummaryProblems <- function(object) {
    e <- S4Vectors::mcols(object@edges)
    c1 <- object@cell_summary_1
    c2 <- object@cell_summary_2
    a <- object@assignment_summary
    statuses <- c("assigned", "unassigned", "ambiguous")
    totals <- c(.countTotal(a$n_molecules), .countTotal(e$weight),
                .countTotal(c1$assigned_mass), .countTotal(c2$assigned_mass))
    if (!all(c("status_1", "status_2", "n_molecules") %in% colnames(a)) ||
        !identical(a$status_1, rep(statuses, each = 3L)) ||
        !identical(a$status_2, rep(statuses, 3L)) ||
        anyNA(totals) || totals[1L] != object@universe_size ||
        totals[2L] != a$n_molecules[1L] ||
        totals[3L] != sum(a$n_molecules[1:3]) ||
        totals[4L] != sum(a$n_molecules[c(1L, 4L, 7L)])) {
        return("Invalid assignment-status conservation summary")
    }
    TRUE
}

.comparisonValidity <- function(object) {
    nm <- object@segmentation_names
    if (!is.character(nm) || length(nm) != 2L || anyNA(nm) ||
        any(!nzchar(nm)) || anyDuplicated(nm)) {
        return("Comparison requires two distinct segmentation names")
    }
    hashes <- c(object@universe_fingerprint, object@assignment_fingerprints)
    if (length(object@universe_fingerprint) != 1L ||
        length(object@assignment_fingerprints) != 2L || anyNA(hashes) ||
        any(!grepl("^[0-9a-f]{64}$", hashes)) ||
        !identical(object@fingerprint_version, "1") ||
        !identical(object@fingerprint_algorithm, "sha256") ||
        !identical(object@schema_version, "1") ||
        length(object@universe_size) != 1L ||
        !.exactCounts(object@universe_size)) {
        return(paste("Invalid comparison fingerprint metadata, schema",
                     "or universe_size"))
    }
    c1 <- object@cell_summary_1
    c2 <- object@cell_summary_2
    for (cells in list(c1, c2)) {
        problems <- .cellSummaryProblems(cells)
        if (!isTRUE(problems)) return(problems)
    }
    h <- object@edges
    problems <- .graphStructureProblems(h, c1, c2)
    if (!isTRUE(problems)) return(problems)
    problems <- .directionalEdgeProblems(h, c1, c2)
    if (!isTRUE(problems)) return(problems)
    problems <- .componentProblems(object@components, c1, c2)
    if (!isTRUE(problems)) return(problems)
    if (any(c1$component_id[S4Vectors::from(h)] !=
            c2$component_id[S4Vectors::to(h)])) {
        return("Edge endpoints have different component IDs")
    }
    .assignmentSummaryProblems(object)
}

#' Complete descriptive correspondence between two segmentations
#'
#' Represent a complete positive shared-molecule graph, independent assigned
#' denominators, directional dominance and structural components. Results do
#' not establish biological truth or rank segmentation accuracy.
#' @rdname SegmentationComparison
#' @md
#' @details Create results with [compareSegmentations()]; there is no public
#'   comparison constructor or replacement API. Query and subject indices in
#'   the Hits returned by correspondenceEdges() refer to the full side-1 and
#'   side-2 census rows. Edges are unique positive pairs sorted by query then
#'   subject. External cell IDs remain lossless in the side tables.
#'
#'   Edge metadata columns are weight, fraction_1, fraction_2, dominant_1,
#'   dominant_2, unique_dominant_1, unique_dominant_2 and reciprocal_dominant.
#'   Directional fractions divide shared weight by the endpoint's complete
#'   assigned mass. Dominance maximizes those fractions. Comparing integer
#'   weights is exactly equivalent because the denominator is constant across
#'   an endpoint's edges, and avoids rounding-induced ties. Every exact tie is
#'   retained. Reciprocity requires a unique maximum on both endpoints.
#'
#'   Cell tables contain cell_id, assigned_mass, matched_mass, unmatched_mass,
#'   partner_count, dominance_count, component_id and matched_fraction, in
#'   original census order. Matched mass sums incident weights; unmatched mass
#'   is assigned minus matched. Zero-molecule and zero-correspondence cells
#'   remain explicit. No edges means zero dominance_count; one maximum means
#'   unique dominance; more than one means a tie. matched_fraction is NA for
#'   zero assigned mass, and otherwise matched_mass/assigned_mass.
#'
#'   Components use the complete positive graph including isolates. The table
#'   contains component_id, n_segmentation_1, n_segmentation_2 and
#'   component_class. Cardinalities (1,1), (1,>1), (>1,1), (>1,>1) and isolated
#'   (1,0)/(0,1) give one_to_one, split, merge, many_to_many and unmatched.
#'   These are structural summaries, not biological split/merge events or
#'   anchor topology. IDs follow first encountered node in side-1 then side-2
#'   census order; census permutations may change component numbering.
#'
#'   Assignment accounting has nine rows with status_1, status_2 and
#'   n_molecules.
#'   Each side uses assigned, unassigned, ambiguous order; side 1 varies
#'   slowest.
#'   Every molecule contributes once. The assigned/assigned count equals total
#'   edge weight and matched mass on either side. Unresolved statuses are not
#'   forced onto edges or conflated with unmatched assigned mass.
#'
#'   Counts are exact whole-valued numeric values at most 2^53. Comparisons
#'   retain internal source bindings, coordinate declarations and nested set/
#'   segmentation provenance with unknown values preserved. They contain no
#'   molecule universe, molecule assignment vectors, geometry or states.
#'   Structural validity does not recompute fingerprints. Standard R copy and
#'   serialization semantics apply, including reference-valued ancillary data.
#' @param x A SegmentationComparison.
#' @param object A SegmentationComparison to display or summarize.
#' @param segmentation Exactly one of the two compared segmentation names.
#' @param ... Additional summary arguments; currently unused.
#' @return correspondenceEdges() returns Hits. cellSummary(), componentSummary()
#'   and assignmentSummary() return DataFrames. segmentationNames() returns
#'   the two names in comparison order; provenance() returns a SimpleList with
#'   set and named segmentations entries. summary() returns a SimpleList with
#'   segmentation_names, universe_size, cell_summary_1, cell_summary_2,
#'   components and assignment_summary. show() displays structural counts.
#' @seealso [compareSegmentations()], [SegmentationSet()]
#' @examples
#' m <- S4Vectors::DataFrame(molecule_id = c("m1", "m2"),
#'                          feature_id = c("g", "g"))
#' p <- SegmentationPartition(m, "left",
#'     S4Vectors::DataFrame(cell_id = c("x", "empty")),
#'     m$molecule_id, c("assigned", "assigned"), c("x", "x"))
#' q <- SegmentationPartition(m, "right",
#'     S4Vectors::DataFrame(cell_id = c("a", "b")),
#'     m$molecule_id, c("assigned", "assigned"), c("a", "b"))
#' cmp <- compareSegmentations(
#'     SegmentationSet(m, list(left = p, right = q), "specimen"),
#'     "left", "right")
#' cmp
#' correspondenceEdges(cmp)
#' cellSummary(cmp, "left")
#' componentSummary(cmp)
#' assignmentSummary(cmp)
#' segmentationNames(cmp)
#' provenance(cmp)
#' summary(cmp)
#' @exportClass SegmentationComparison
methods::setClass("SegmentationComparison", slots = c(
    segmentation_names = "character", cell_summary_1 = "DataFrame",
    cell_summary_2 = "DataFrame", edges = "Hits", components = "DataFrame",
    assignment_summary = "DataFrame", universe_size = "numeric",
    universe_fingerprint = "character", assignment_fingerprints = "character",
    fingerprint_version = "character", fingerprint_algorithm = "character",
    coordinate_info = "SimpleList", provenance = "SimpleList",
    schema_version = "character", metadata = "SimpleList"),
    validity = .comparisonValidity)

.newComparison <- function(...) methods::new("SegmentationComparison", ...)

#' @rdname SegmentationComparison
#' @export
methods::setMethod("show", "SegmentationComparison", function(object) {
    cat("SegmentationComparison:", paste(object@segmentation_names,
                                         collapse = " / "), "\n")
    cat(" ", object@universe_size, "molecules;", nrow(object@cell_summary_1),
        "/", nrow(object@cell_summary_2), "cells;", length(object@edges),
        "positive edges;", nrow(object@components), "components\n")
})

#' @rdname SegmentationComparison
#' @export
methods::setMethod("summary", "SegmentationComparison", function(object, ...) {
    S4Vectors::SimpleList(segmentation_names = object@segmentation_names,
        universe_size = object@universe_size,
        cell_summary_1 = object@cell_summary_1,
        cell_summary_2 = object@cell_summary_2, components = object@components,
        assignment_summary = object@assignment_summary)
})
