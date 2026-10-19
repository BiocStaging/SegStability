.supportRequire <- function(report) {
    blocked <- report$status %in% c("FAIL", "UNKNOWN")
    if (any(blocked)) {
        stop("Support-analysis prerequisite did not PASS: ",
            paste(paste0(report$check[blocked], " [", report$status[blocked],
                         "]: ", report$message[blocked]), collapse = "; "),
            call. = FALSE)
    }
    report
}

.supportInputs <- function(thresholds, direction, names) {
    threshold_check <- .validationCheck("support.thresholds", {
        if (!.baseMeasurement(thresholds) || !is.null(dim(thresholds)) ||
            !length(thresholds) || anyNA(thresholds) ||
            any(!is.finite(thresholds) | thresholds < 0 | thresholds > 1) ||
            anyDuplicated(thresholds)) {
            stop("thresholds must be unique finite plain fractions in [0,1]")
        }
        TRUE
    })
    .supportRequire(threshold_check)
    direction_check <- .validationCheck("support.direction", {
        if (!is.null(direction) &&
            (!identical(class(direction), "character") ||
             !is.null(dim(direction)) || length(direction) != 1L ||
             anyNA(direction) || !direction %in% names)) {
            stop("direction must be NULL or one exact segmentation name")
        }
        TRUE
    })
    .supportRequire(direction_check)
    list(thresholds = as.numeric(thresholds),
        sides = if (is.null(direction)) 1:2 else match(direction, names),
        report = .validationBind(list(threshold_check, direction_check)))
}

.supportStructure <- function(original, index, weight, keep) {
    n <- nrow(original)
    retained <- .sumByIndex(index[keep], weight[keep], n)
    count <- tabulate(index[keep], n)
    list(retained_mass = retained,
        removed_mass = original$matched_mass - retained,
        edge_count = count, removed_edge_count = original$partner_count - count,
        retained_mass_fraction = .stateRatio(retained, original$matched_mass),
        retained_assigned_fraction =
            .stateRatio(retained, original$assigned_mass))
}

.supportStates <- function(index, weight, keep, resolved, same, values, cells) {
    n <- nrow(cells)
    fields <- c("resolved_mass", "same_state_mass", "different_state_mass",
        "unresolved_label_mass", "conditional_stability",
        "same_state_fraction_of_matched_mass", "total_mass_stability",
        "endpoint_unavailable_cells", "total_estimable_cells",
        "total_nonestimable_cells", "conditional_estimable_cells",
        "conditional_nonestimable_cells", "mean_cell_total_mass_stability",
        "mean_cell_conditional_stability")
    if (is.null(values)) {
        return(stats::setNames(rep(list(rep(NA_real_, n)), length(fields)),
                               fields))
    }
    r <- .sumByIndex(index[keep & resolved], weight[keep & resolved], n)
    s <- .sumByIndex(index[keep & same], weight[keep & same], n)
    m <- cells$selected_retained_mass
    total <- .stateRatio(s, cells$assigned_mass)
    conditional <- .stateRatio(s, r)
    list(resolved_mass = r, same_state_mass = s, different_state_mass = r - s,
        unresolved_label_mass = m - r, conditional_stability = conditional,
        same_state_fraction_of_matched_mass = .stateRatio(s, m),
        total_mass_stability = total,
        endpoint_unavailable_cells = as.integer(is.na(as.character(values))),
        total_estimable_cells = as.integer(cells$assigned_mass > 0),
        total_nonestimable_cells = as.integer(cells$assigned_mass == 0),
        conditional_estimable_cells = as.integer(r > 0),
        conditional_nonestimable_cells = as.integer(r == 0),
        mean_cell_total_mass_stability = total,
        mean_cell_conditional_stability = conditional)
}

.supportCellTables <- function(comparison, weight, keep, labels) {
    indices <- list(S4Vectors::from(comparison@edges),
                    S4Vectors::to(comparison@edges))
    resolved <- same <- rep(FALSE, length(weight))
    if (!is.null(labels)) {
        a <- as.character(labels[[1L]])[indices[[1L]]]
        b <- as.character(labels[[2L]])[indices[[2L]]]
        resolved <- !is.na(a) & !is.na(b)
        same[resolved] <- a[resolved] == b[resolved]
    }
    tables <- list(comparison@cell_summary_1, comparison@cell_summary_2)
    for (side in 1:2) {
        c <- tables[[side]]
        structural <- .supportStructure(c, indices[[side]], weight, keep)
        for (field in names(structural)) {
            c[[paste0("selected_", field)]] <- structural[[field]]
        }
        values <- if (is.null(labels)) NULL else labels[[side]]
        state <- .supportStates(indices[[side]], weight, keep, resolved, same,
                                values, c)
        for (field in names(state)) {
            c[[paste0("selected_", field)]] <- state[[field]]
        }
        if (!is.null(labels)) c$state <- values
        tables[[side]] <- c
    }
    out <- S4Vectors::SimpleList(tables)
    names(out) <- comparison@segmentation_names
    out
}

.supportSummary <- function(cells, states_supplied) {
    m <- .countTotal(cells$selected_retained_mass)
    assigned <- .countTotal(cells$assigned_mass)
    out <- S4Vectors::DataFrame(
        original_edges = as.integer(sum(cells$partner_count)),
        retained_edges = as.integer(sum(cells$selected_edge_count)),
        removed_edges = as.integer(sum(cells$selected_removed_edge_count)),
        original_matched_mass = .countTotal(cells$matched_mass),
        retained_mass = m,
        excluded_matched_mass = .countTotal(cells$selected_removed_mass),
        retained_mass_fraction =
            .stateRatio(m, .countTotal(cells$matched_mass)),
        assigned_mass = assigned,
        unmatched_mass = .countTotal(cells$unmatched_mass),
        retained_assigned_fraction = .stateRatio(m, assigned),
        n_cells = nrow(cells),
        zero_molecule_cells = sum(cells$assigned_mass == 0),
        zero_correspondence_cells = sum(cells$matched_mass == 0),
        retained_cells = sum(cells$selected_retained_mass > 0),
        selected_zero_correspondence_cells =
            sum(cells$selected_retained_mass == 0))
    fields <- c("resolved_mass", "same_state_mass", "different_state_mass",
                "unresolved_label_mass", "endpoint_unavailable_cells",
                "total_estimable_cells", "total_nonestimable_cells",
                "conditional_estimable_cells", "conditional_nonestimable_cells")
    for (field in fields) {
        out[[field]] <- if (states_supplied) {
            .countTotal(cells[[paste0("selected_", field)]])
        } else NA_real_
    }
    out$conditional_stability <- .stateRatio(out$same_state_mass,
                                            out$resolved_mass)
    out$same_state_fraction_of_matched_mass <-
        .stateRatio(out$same_state_mass, m)
    out$total_mass_stability <- .stateRatio(out$same_state_mass, assigned)
    for (kind in c("total_mass", "conditional")) {
        estimable <- if (kind == "total_mass") cells$assigned_mass > 0 else
            !is.na(cells$selected_resolved_mass) &
                cells$selected_resolved_mass > 0
        out[[paste0("mean_cell_", kind, "_stability")]] <-
            if (states_supplied && any(estimable)) {
                field <- paste0("selected_", kind, "_stability")
                mean(cells[[field]][estimable])
            } else NA_real_
    }
    out
}

.supportComponentSummary <- function(cells, groups, states_supplied) {
    columns <- as.list(cells)
    n <- length(groups)
    mass <- function(field) {
        values <- columns[[field]]
        unname(vapply(groups, function(i) .countTotal(values[i]), numeric(1)))
    }
    count <- function(values) {
        unname(vapply(groups, function(i) as.integer(sum(values[i])), integer(1)))
    }
    m <- mass("selected_retained_mass")
    assigned <- mass("assigned_mass")
    matched <- mass("matched_mass")
    out <- list(original_edges = count(columns$partner_count),
        retained_edges = count(columns$selected_edge_count),
        removed_edges = count(columns$selected_removed_edge_count),
        original_matched_mass = matched, retained_mass = m,
        excluded_matched_mass = mass("selected_removed_mass"),
        retained_mass_fraction = .stateRatio(m, matched),
        assigned_mass = assigned, unmatched_mass = mass("unmatched_mass"),
        retained_assigned_fraction = .stateRatio(m, assigned),
        n_cells = unname(lengths(groups)),
        zero_molecule_cells = count(columns$assigned_mass == 0),
        zero_correspondence_cells = count(columns$matched_mass == 0),
        retained_cells = count(columns$selected_retained_mass > 0),
        selected_zero_correspondence_cells =
            count(columns$selected_retained_mass == 0))
    fields <- c("resolved_mass", "same_state_mass", "different_state_mass",
        "unresolved_label_mass", "endpoint_unavailable_cells",
        "total_estimable_cells", "total_nonestimable_cells",
        "conditional_estimable_cells", "conditional_nonestimable_cells")
    for (field in fields) {
        out[[field]] <- if (states_supplied) mass(paste0("selected_", field)) else
            rep(NA_real_, n)
    }
    out$conditional_stability <- if (states_supplied) {
        .stateRatio(out$same_state_mass, out$resolved_mass)
    } else rep(NA_real_, n)
    out$same_state_fraction_of_matched_mass <-
        .stateRatio(out$same_state_mass, m)
    out$total_mass_stability <- .stateRatio(out$same_state_mass, assigned)
    for (kind in c("total_mass", "conditional")) {
        eligible <- if (kind == "total_mass") columns$assigned_mass > 0 else
            !is.na(columns$selected_resolved_mass) &
                columns$selected_resolved_mass > 0
        values <- columns[[paste0("selected_", kind, "_stability")]]
        out[[paste0("mean_cell_", kind, "_stability")]] <-
            if (states_supplied) {
                unname(vapply(groups, function(i) {
                    selected <- i[eligible[i]]
                    if (length(selected)) mean(values[selected]) else NA_real_
                }, numeric(1)))
            } else rep(NA_real_, n)
    }
    out
}

.supportComponents <- function(comparison, cells, states_supplied) {
    out <- comparison@components
    n <- nrow(out)
    for (side in 1:2) {
        c <- cells[[side]]
        groups <- split(seq_len(nrow(c)), factor(c$component_id,
                                                levels = seq_len(n)))
        # Reduce aligned vectors within the original groups, including groups
        # with no cells on this side, before materializing any component columns.
        table <- .supportComponentSummary(c, groups, states_supplied)
        for (field in names(table)) {
            out[[paste0("selected_", field, "_", side)]] <- table[[field]]
        }
        out[[paste0("assigned_mass_", side)]] <- table$assigned_mass
        out[[paste0("unmatched_mass_", side)]] <- table$unmatched_mass
        if (side == 1L) {
            # Every edge belongs to one original component on both endpoints.
            fields <- c("original_edges", "original_matched_mass",
                "retained_edges", "removed_edges", "retained_mass",
                "excluded_matched_mass", "retained_mass_fraction",
                "resolved_mass", "same_state_mass", "different_state_mass",
                "unresolved_label_mass", "conditional_stability",
                "same_state_fraction_of_matched_mass")
            for (field in fields) {
                key <- switch(field, original_edges = "original_edges",
                    original_matched_mass = "original_matched_mass",
                    retained_edges = "selected_edge_count",
                    removed_edges = "selected_removed_edge_count",
                    excluded_matched_mass = "selected_removed_mass",
                    paste0("selected_", field))
                out[[key]] <- table[[field]]
            }
        }
    }
    out
}

#' Describe directional support sensitivity of complete correspondence
#'
#' Measure support retention and optional categorical agreement using the
#' immutable complete comparison and its original independent denominators.
#' @param comparison A SegmentationComparison from compareSegmentations().
#' @param thresholds Nonempty plain integer or numeric vector of unique finite
#'   fractions in [0,1]. Supplied order is retained.
#' @param direction NULL for both segmentation names in comparison order, or
#'   one exact segmentation name.
#' @param states NULL, or a named list following the stateStability() mapping
#'   contract. States do not affect structural support selection or accounting.
#' @details
#' For each direction select edges with weight/original assigned_mass greater
#' than or equal to the threshold. The directions are separate estimands, not
#' an AND/OR selection. No molecule rematching or stateStability() call occurs.
#' Original graph metadata, dominance/ties, reciprocity, component identities,
#' assignment accounting and cell censuses remain unchanged. Selected support
#' is summarized within original components; no component reclassification,
#' geometry selection or anchor topology is performed. Selected edges are
#' internal views, never a new primary comparison.
#'
#' Retained mass fraction uses all original matched mass. Retained assigned
#' fraction and total_mass_stability use all originally assigned mass on the
#' selected side. Original unmatched mass excludes only original nonmatches;
#' excluded_matched_mass reports the separate cost of support filtering.
#' Supplied states yield selected resolved mass R, same-state mass S,
#' different-state mass R-S and unresolved-label mass retained_mass-R.
#' Conditional stability is S/R; total-mass stability is S/N_original.
#' Missing labels are unresolved, never disagreements. With states=NULL all
#' state-dependent fields are NA. Zero denominators yield NA; positive
#' denominators with zero numerator yield zero. Equal-cell total means use
#' original assigned_mass>0; conditional means use selected resolved_mass>0.
#'
#' Required support checks must PASS. Unavailable source-level fingerprint
#' recomputation, external census completeness and authenticity evidence stays
#' UNKNOWN and does not block otherwise valid support analysis. Validation
#' metadata does not certify source authenticity or independently reverified
#' source integrity. These results are descriptive, not accuracy rankings.
#' @return A DataFrame, with direction-major rows then supplied threshold order.
#' Scalar fields are direction, threshold, states_supplied, original_edges,
#' retained_edges, removed_edges, original_matched_mass, retained_mass,
#' excluded_matched_mass, retained_mass_fraction, assigned_mass, unmatched_mass,
#' retained_assigned_fraction, n_cells, zero_molecule_cells,
#' zero_correspondence_cells (original), retained_cells and
#' selected_zero_correspondence_cells. State fields are resolved_mass,
#' same_state_mass, different_state_mass, unresolved_label_mass,
#' conditional_stability, same_state_fraction_of_matched_mass,
#' total_mass_stability, endpoint_unavailable_cells, total_estimable_cells,
#' total_nonestimable_cells, conditional_estimable_cells,
#' conditional_nonestimable_cells, mean_cell_total_mass_stability and
#' mean_cell_conditional_stability.
#'
#' The cells list column holds a SimpleList keyed by both segmentation names,
#' retaining original census fields and appending selected_retained_mass,
#' selected_removed_mass, selected_edge_count, selected_removed_edge_count,
#' selected_retained_mass_fraction, selected_retained_assigned_fraction and
#' selected_ state fields above (counts become per-cell indicators; means become
#' per-cell ratios). Supplied state vectors retain their types in state.
#' The components list column retains original component rows/classes, with
#' original_edges, original_matched_mass, pooled selected support/state fields
#' and directional scalar summary fields prefixed selected_ and suffixed _1/_2.
#' Pooled selected counts are selected_edge_count/selected_removed_edge_count;
#' pooled removed mass is selected_removed_mass. Original assigned_mass and
#' unmatched_mass are also included with _1/_2 suffixes. Both nested tables
#' preserve zero-support census members. metadata(result)$validation holds
#' DataFrame(check,status,severity,message); support. checks are required,
#' source/provenance diagnostics are contextual and can remain UNKNOWN.
#' @seealso compareSegmentations(), stateStability(), componentSummary()
#' @examples
#' m <- S4Vectors::DataFrame(molecule_id = c("m1", "m2"),
#'                          feature_id = c("g", "g"))
#' p <- SegmentationPartition(m, "left",
#'     S4Vectors::DataFrame(cell_id = "x"), m$molecule_id,
#'     c("assigned", "assigned"), c("x", "x"))
#' q <- SegmentationPartition(m, "right",
#'     S4Vectors::DataFrame(cell_id = c("a", "b")), m$molecule_id,
#'     c("assigned", "assigned"), c("a", "b"))
#' cmp <- compareSegmentations(
#'     SegmentationSet(m, list(left = p, right = q), "specimen"),
#'     "left", "right")
#' result <- supportSensitivity(cmp, thresholds = c(0, .5, 1))
#' result[, c("direction", "threshold", "retained_mass")]
#' result$cells[[1]]$left
#' result$components[[1]]
#' S4Vectors::metadata(result)$validation
#' @export
supportSensitivity <- function(comparison,
    thresholds = c(0, 0.01, 0.025, 0.05, 0.10, 0.20),
    direction = NULL, states = NULL) {
    prerequisites <- .stateComparisonChecks(comparison, "support",
                                             .supportRequire)
    input <- .supportInputs(thresholds, direction,
                            comparison@segmentation_names)
    labels <- NULL
    reports <- list(prerequisites, input$report)
    if (!is.null(states)) {
        mapping <- .stateMappings(states, comparison)
        labels <- mapping$values
        mapping$report$check <- sub("^state\\.", "support.states.",
                                    mapping$report$check)
        reports[[3L]] <- mapping$report
    }
    e <- S4Vectors::mcols(comparison@edges)
    rows <- vector("list", length(input$sides) * length(input$thresholds))
    k <- 0L
    for (side in input$sides) {
        for (threshold in input$thresholds) {
            keep <- e[[paste0("fraction_", side)]] >= threshold
            cells <- .supportCellTables(comparison, e$weight, keep, labels)
            components <- .supportComponents(comparison, cells,
                                              !is.null(states))
            row <- S4Vectors::DataFrame(
                direction = comparison@segmentation_names[side],
                threshold = threshold, states_supplied = !is.null(states))
            row <- cbind(row, .supportSummary(cells[[side]], !is.null(states)))
            row$cells <- I(list(cells))
            row$components <- I(list(components))
            k <- k + 1L
            rows[[k]] <- row
        }
    }
    out <- do.call(rbind, rows)
    reports[[length(reports) + 1L]] <- .stateSourceEvidence(comparison)
    S4Vectors::metadata(out)$validation <- .validationBind(reports)
    out
}
