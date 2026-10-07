.stateRequire <- function(report) {
    failed <- report$status != "PASS"
    if (any(failed)) {
        stop("State-analysis prerequisite did not PASS: ",
            paste(paste0(report$check[failed], " [", report$status[failed],
                         "]: ", report$message[failed]), collapse = "; "),
            call. = FALSE)
    }
    report
}

.stateComparisonChecks <- function(comparison) {
    input <- .validationCheck("state.comparison.input",
        if (methods::is(comparison, "SegmentationComparison")) TRUE else
            "Expected a SegmentationComparison")
    .stateRequire(input)
    columns <- .validationCheck("state.comparison.columns", {
        tables <- list(comparison@cell_summary_1, comparison@cell_summary_2,
            comparison@components, comparison@assignment_summary,
            S4Vectors::mcols(comparison@edges))
        if (any(vapply(tables, function(x) anyDuplicated(colnames(x)) > 0L,
                       logical(1)))) {
            stop("Scientific comparison tables require unique column names")
        }
        TRUE
    })
    .stateRequire(columns)
    structure <- .validationCheck("state.comparison.structure",
        methods::validObject(comparison, test = TRUE, complete = TRUE))
    .stateRequire(structure)
    graph <- .validationCheck("state.comparison.complete_graph", {
        h <- comparison@edges
        q <- S4Vectors::from(h)
        s <- S4Vectors::to(h)
        e <- S4Vectors::mcols(h)
        c1 <- comparison@cell_summary_1
        c2 <- comparison@cell_summary_2
        for (side in 1:2) {
            idx <- if (side == 1L) q else s
            cells <- if (side == 1L) c1 else c2
            d <- .directionalDominance(idx, e$weight, nrow(cells))
            if (!identical(d$dominant, e[[paste0("dominant_", side)]]) ||
                !identical(d$unique,
                           e[[paste0("unique_dominant_", side)]]) ||
                any(d$count != cells$dominance_count)) {
                stop("Complete-graph dominance or exact ties are inconsistent")
            }
        }
        components <- .graphComponents(q, s, nrow(c1), nrow(c2))
        if (!identical(components$side1, c1$component_id) ||
            !identical(components$side2, c2$component_id)) {
            stop("Stored components differ from complete connected components")
        }
        TRUE
    })
    .stateRequire(graph)
    .validationBind(list(input, columns, structure, graph))
}

.stateValues <- function(x) {
    if (!is.null(dim(x)) ||
        !(identical(class(x), "character") ||
          identical(class(x), "integer") ||
          identical(class(x), "factor") ||
          identical(class(x), c("ordered", "factor")))) {
        stop("state must be a character, factor or integer categorical vector",
             call. = FALSE)
    }
    x
}

.stateAlign <- function(mapping, cells) {
    if (methods::is(mapping, "DataFrame")) {
        if (anyDuplicated(colnames(mapping)) ||
            !all(c("cell_id", "state") %in% colnames(mapping))) {
            stop("State DataFrame requires unique cell_id and state columns")
        }
        ids <- mapping$cell_id
        values <- .stateValues(mapping$state)
    } else {
        values <- .stateValues(mapping)
        if (.idFamily(cells$cell_id) != "character") {
            stop("Non-character cell IDs require a typed state DataFrame")
        }
        ids <- names(mapping)
        if (is.null(ids)) stop("State vector requires cell-ID names")
    }
    .validateIds(ids, "state cell_id")
    if (length(values) != length(ids)) stop("State mapping dimensions differ")
    index <- .matchIds(cells$cell_id, ids)
    if (anyNA(index) || anyNA(.matchIds(ids, cells$cell_id))) {
        stop("State mapping must cover exactly the full cell census")
    }
    # Names encode alignment, not categorical values; retain storage and levels.
    unname(values[index])
}

.stateMappings <- function(states, comparison) {
    nm <- comparison@segmentation_names
    shape <- .validationCheck("state.mapping.names", {
        keys <- names(states)
        if (!is.list(states) || is.object(states) || length(states) != 2L ||
            is.null(keys) || anyNA(keys) || anyDuplicated(keys) ||
            !setequal(keys, nm)) {
            stop("states must be a named list keyed by both segmentation names")
        }
        TRUE
    })
    .stateRequire(shape)
    aligned <- vector("list", 2L)
    reports <- list(shape)
    for (i in 1:2) {
        cells <- if (i == 1L) comparison@cell_summary_1 else
            comparison@cell_summary_2
        reports[[i + 1L]] <- .validationCheck(
            paste0("state.mapping[", i, "].alignment"), {
                aligned[[i]] <- .stateAlign(states[[nm[i]]], cells)
                TRUE
            })
        .stateRequire(reports[[i + 1L]])
    }
    list(values = aligned, report = .validationBind(reports))
}

.stateRatio <- function(numerator, denominator) {
    result <- rep(NA_real_, length(denominator))
    positive <- denominator > 0
    result[positive] <- numerator[positive] / denominator[positive]
    result
}

.stateMassSummary <- function(matched, resolved, same) {
    S4Vectors::DataFrame(matched_mass = matched, resolved_mass = resolved,
        same_state_mass = same, different_state_mass = resolved - same,
        unresolved_label_mass = matched - resolved,
        conditional_stability = .stateRatio(same, resolved),
        same_state_fraction_of_matched_mass = .stateRatio(same, matched))
}

.stateCells <- function(cells, state, index, weight, resolved, same) {
    cells$state <- state
    cells$resolved_mass <- .sumByIndex(index[resolved], weight[resolved],
                                      nrow(cells))
    cells$same_state_mass <- .sumByIndex(index[same], weight[same], nrow(cells))
    cells$different_state_mass <- cells$resolved_mass - cells$same_state_mass
    cells$unresolved_label_mass <- cells$matched_mass - cells$resolved_mass
    cells$total_mass_stability <- .stateRatio(cells$same_state_mass,
                                            cells$assigned_mass)
    cells$conditional_stability <- .stateRatio(cells$same_state_mass,
                                             cells$resolved_mass)
    cells
}

.stateCellSummary <- function(cells) {
    out <- .stateMassSummary(.countTotal(cells$matched_mass),
        .countTotal(cells$resolved_mass), .countTotal(cells$same_state_mass))
    out$assigned_mass <- .countTotal(cells$assigned_mass)
    out$unmatched_mass <- out$assigned_mass - out$matched_mass
    out$total_mass_stability <- .stateRatio(out$same_state_mass,
                                          out$assigned_mass)
    total <- cells$assigned_mass > 0
    conditional <- cells$resolved_mass > 0
    out$n_cells <- nrow(cells)
    out$zero_molecule_cells <- sum(!total)
    out$zero_correspondence_cells <- sum(cells$matched_mass == 0)
    # An NA factor level can have a nonmissing code; availability uses labels.
    out$endpoint_unavailable_cells <- sum(is.na(as.character(cells$state)))
    out$total_estimable_cells <- sum(total)
    out$total_nonestimable_cells <- sum(!total)
    out$conditional_estimable_cells <- sum(conditional)
    out$conditional_nonestimable_cells <- sum(!conditional)
    out$mean_cell_total_mass_stability <- if (any(total)) {
        mean(cells$total_mass_stability[total])
    } else NA_real_
    out$mean_cell_conditional_stability <- if (any(conditional)) {
        mean(cells$conditional_stability[conditional])
    } else NA_real_
    out
}

.stateComponents <- function(comparison, cells, edge_mass) {
    out <- comparison@components
    n <- nrow(out)
    edge_component <- cells[[1L]]$component_id[
        S4Vectors::from(comparison@edges)]
    mass <- lapply(edge_mass, function(x) .sumByIndex(edge_component, x, n))
    pooled <- .stateMassSummary(mass[[1L]], mass[[2L]], mass[[3L]])
    for (field in colnames(pooled)) out[[field]] <- pooled[[field]]
    for (side in 1:2) {
        c <- cells[[side]]
        # Group positions once: isolates and empty opposite sides stay explicit.
        groups <- split(seq_len(nrow(c)), factor(c$component_id,
                                                levels = seq_len(n)))
        rows <- lapply(groups, function(i) .stateCellSummary(c[i, ]))
        summary <- if (length(rows)) do.call(rbind, rows) else
            .stateCellSummary(c[integer(), ])[integer(), ]
        fields <- setdiff(colnames(summary), colnames(pooled))
        for (field in fields) {
            out[[paste0(field, "_", side)]] <- summary[[field]]
        }
    }
    out
}

.stateGlobal <- function(comparison, directional) {
    d <- directional
    out <- .stateMassSummary(d$matched_mass[1L], d$resolved_mass[1L],
                             d$same_state_mass[1L])
    for (i in 1:2) {
        out[[paste0("segmentation_", i)]] <- comparison@segmentation_names[i]
        for (field in c("assigned_mass", "unmatched_mass",
                        "total_mass_stability")) {
            out[[paste0(field, "_", i)]] <- d[[field]][i]
        }
    }
    a <- comparison@assignment_summary
    ambiguous1 <- a$status_1 == "ambiguous"
    ambiguous2 <- a$status_2 == "ambiguous"
    out$ambiguous_mass_1 <- .countTotal(a$n_molecules[ambiguous1])
    out$ambiguous_mass_2 <- .countTotal(a$n_molecules[ambiguous2])
    # Each contingency row is disjoint molecule mass, including both-ambiguous.
    out$ambiguous_mass_union <- .countTotal(
        a$n_molecules[ambiguous1 | ambiguous2])
    out
}

.stateSourceEvidence <- function(comparison) {
    .validationBind(list(
        .validationUnknown("source.fingerprint_recomputation",
            "Source molecules/assignments are not retained; not recomputed"),
        .validationUnknown("source.census_completeness",
            "External source census completeness cannot be established"),
        validateProvenance(comparison)))
}

#' Describe categorical state stability across complete correspondence
#'
#' Quantify same-state, different-state and unresolved shared-molecule mass
#' without ranking segmentation accuracy or claiming biological truth.
#' @param comparison A SegmentationComparison from compareSegmentations().
#' @param states A base named list keyed by exactly both compared segmentation
#'   names. Each element is a named character, factor or integer categorical
#'   vector, or an S4Vectors DataFrame with cell_id and state columns. Each
#'   mapping must cover its full census exactly. Additional DataFrame columns
#'   are not used. Named vectors require character census IDs; other ID
#'   families require a DataFrame with the same lossless cell-ID family.
#' @details
#' Missing mapping cells, unknown or duplicate IDs and unsupported state types
#' are errors. A present cell with typed NA state or an NA factor label is
#' unresolved. Factors compare by labels; integer labels compare by exact
#' decimal text internally, so 1L
#' equals "1" but differs from "01". Original values, types and factor levels
#' are retained in public cell and edge tables. No normalization or synonym
#' mapping is performed; an empty character label is an available category.
#'
#' Let M be all edge weight, R weight with both endpoint labels available, and
#' S resolved same-state weight. Different-state mass is R-S and unresolved
#' label mass is M-R. Missing labels are neither agreements nor disagreements.
#' Conditional stability is S/R; the same-state matched fraction is S/M.
#' Directional total-mass stability is S/N using each complete independently
#' assigned population, including assignments unresolved on the other side.
#' Each cell uses all its incident edges. Zero denominators return NA; a zero
#' numerator with positive denominator returns zero. Equal-cell total means
#' use N_i>0; conditional means use R_i>0, with population counts reported.
#'
#' State-analysis prerequisites are explicitly checked and must PASS, including
#' structure, complete-graph accounting, dominance, components and state
#' alignment. Unavailable source fingerprint recomputation, external census
#' completeness and authenticity remain UNKNOWN and do not block an otherwise
#' valid state calculation. Returned diagnostics do not attest independently
#' reverified source integrity. Provenance diagnostics are contextual evidence,
#' separate from required state checks. No source data or validation cache is
#' inferred, and no molecule rematching or support/topology selection occurs.
#' Ambiguity is reported from the molecule-status contingency: side-specific
#' counts and the distinct-molecule union, counting both-ambiguous molecules
#' once. Ambiguous assignments never acquire forced correspondence edges.
#' @return A SimpleList with elements in this order:
#' \describe{
#' \item{global}{One-row DataFrame with matched_mass (M), resolved_mass (R),
#' same_state_mass (S), different_state_mass, unresolved_label_mass,
#' conditional_stability, same_state_fraction_of_matched_mass,
#' segmentation_1/2, assigned_mass_1/2, unmatched_mass_1/2,
#' total_mass_stability_1/2, ambiguous_mass_1/2 and ambiguous_mass_union.}
#' \item{directional}{Two-row DataFrame in comparison order, with segmentation,
#' the unsuffixed masses/ratios above except ambiguity, n_cells,
#' zero_molecule_cells, zero_correspondence_cells, endpoint_unavailable_cells,
#' total_estimable_cells, total_nonestimable_cells,
#' conditional_estimable_cells, conditional_nonestimable_cells,
#' mean_cell_total_mass_stability and mean_cell_conditional_stability.}
#' \item{cells}{SimpleList keyed by segmentation name, retaining original census
#' tables/order and appending state, resolved_mass, same_state_mass,
#' different_state_mass, unresolved_label_mass, total_mass_stability and
#' conditional_stability.}
#' \item{edges}{Hits retaining original indices and metadata, appending state_1,
#' state_2, component_id, resolved, same_state (NA if unresolved),
#' resolved_mass,
#' same_state_mass, different_state_mass and unresolved_label_mass.}
#' \item{components}{DataFrame retaining original component rows/classes and
#' appending pooled M/R/S masses and ratios. Assigned/unmatched mass,
#' total-mass stability and directional cell-count/equal-cell-mean fields are
#' appended with _1/_2 suffixes. No biological component class is inferred.}
#' \item{validation}{DataFrame(check,status,severity,message). Required checks
#' have state. prefixes and PASS; unavailable evidence remains UNKNOWN.}
#' }
#' @seealso compareSegmentations(), cellSummary(), assignmentSummary(),
#'   validateProvenance()
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
#' result <- stateStability(cmp,
#'     list(left = c(x = "s"), right = c(a = "s", b = NA_character_)))
#' result$global
#' result$cells$left
#' result$validation
#' @export
stateStability <- function(comparison, states) {
    prerequisites <- .stateComparisonChecks(comparison)
    mapping <- .stateMappings(states, comparison)
    h <- comparison@edges
    q <- S4Vectors::from(h)
    s <- S4Vectors::to(h)
    w <- S4Vectors::mcols(h)$weight
    values <- mapping$values
    labels <- lapply(values, as.character)
    resolved <- !is.na(labels[[1L]][q]) & !is.na(labels[[2L]][s])
    same_state <- rep(NA, length(w))
    same_state[resolved] <- labels[[1L]][q[resolved]] ==
        labels[[2L]][s[resolved]]
    same <- resolved & !is.na(same_state) & same_state
    cells <- S4Vectors::SimpleList(
        .stateCells(comparison@cell_summary_1, values[[1L]], q, w,
                    resolved, same),
        .stateCells(comparison@cell_summary_2, values[[2L]], s, w,
                    resolved, same))
    names(cells) <- comparison@segmentation_names
    directional <- rbind(.stateCellSummary(cells[[1L]]),
                         .stateCellSummary(cells[[2L]]))
    directional$segmentation <- comparison@segmentation_names
    e <- S4Vectors::mcols(h)
    e$state_1 <- values[[1L]][q]
    e$state_2 <- values[[2L]][s]
    e$component_id <- cells[[1L]]$component_id[q]
    e$resolved <- resolved
    e$same_state <- same_state
    e$resolved_mass <- w * resolved
    e$same_state_mass <- w * same
    e$different_state_mass <- e$resolved_mass - e$same_state_mass
    e$unresolved_label_mass <- w - e$resolved_mass
    S4Vectors::mcols(h) <- e
    S4Vectors::SimpleList(global = .stateGlobal(comparison, directional),
        directional = directional, cells = cells, edges = h,
        components = .stateComponents(comparison, cells,
            list(w, e$resolved_mass, e$same_state_mass)),
        validation = .validationBind(list(prerequisites, mapping$report,
                                         .stateSourceEvidence(comparison))))
}
