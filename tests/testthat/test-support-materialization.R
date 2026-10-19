# Assisted-by: OpenAI Codex
# Independent small-fixture accounting: ordinary tables, no production reducers.
corrective_support_ratio <- function(a, b) {
    as.numeric(ifelse(b > 0, a / b, NA_real_))
}

corrective_support_summary <- function(cells, has_states) {
    mass <- function(name) as.numeric(sum(cells[[name]]))
    ratio <- corrective_support_ratio
    m <- mass("selected_retained_mass")
    original <- mass("matched_mass")
    own <- mass("assigned_mass")
    out <- list(original_edges = as.integer(sum(cells$partner_count)),
        retained_edges = as.integer(sum(cells$selected_edge_count)),
        removed_edges = as.integer(sum(cells$selected_removed_edge_count)),
        original_matched_mass = original, retained_mass = m,
        excluded_matched_mass = original - m,
        retained_mass_fraction = ratio(m, original), assigned_mass = own,
        unmatched_mass = mass("unmatched_mass"),
        retained_assigned_fraction = ratio(m, own), n_cells = nrow(cells),
        zero_molecule_cells = sum(cells$assigned_mass == 0),
        zero_correspondence_cells = sum(cells$matched_mass == 0),
        retained_cells = sum(cells$selected_retained_mass > 0),
        selected_zero_correspondence_cells = sum(cells$selected_retained_mass == 0))
    for (field in c("resolved_mass", "same_state_mass", "different_state_mass",
                   "unresolved_label_mass", "endpoint_unavailable_cells",
                   "total_estimable_cells", "total_nonestimable_cells",
                   "conditional_estimable_cells", "conditional_nonestimable_cells")) {
        out[[field]] <- if (has_states) mass(paste0("selected_", field)) else NA_real_
    }
    out$conditional_stability <- ratio(out$same_state_mass, out$resolved_mass)
    out$same_state_fraction_of_matched_mass <- ratio(out$same_state_mass, m)
    out$total_mass_stability <- ratio(out$same_state_mass, own)
    for (kind in c("total_mass", "conditional")) {
        eligible <- if (kind == "total_mass") cells$assigned_mass > 0 else
            !is.na(cells$selected_resolved_mass) & cells$selected_resolved_mass > 0
        out[[paste0("mean_cell_", kind, "_stability")]] <-
            if (has_states && any(eligible)) {
                mean(cells[[paste0("selected_", kind, "_stability")]][eligible])
            } else NA_real_
    }
    out
}

corrective_support_cells <- function(w, selected, own, labels, opposite,
                                     groups, ids, has_states) {
    n <- nrow(w)
    ratio <- corrective_support_ratio
    matched <- as.numeric(rowSums(w))
    retained <- as.numeric(rowSums(selected))
    degree <- as.integer(rowSums(w > 0))
    retained_degree <- as.integer(rowSums(selected > 0))
    dominant <- vapply(seq_len(n), function(i) {
        if (any(w[i, ] > 0)) sum(w[i, ] == max(w[i, ])) else 0L
    }, integer(1))
    c <- data.frame(cell_id = ids, assigned_mass = own, matched_mass = matched,
        unmatched_mass = own - matched, partner_count = degree,
        dominance_count = dominant, component_id = groups,
        matched_fraction = ratio(matched, own), stringsAsFactors = FALSE)
    c$selected_retained_mass <- retained
    c$selected_removed_mass <- matched - retained
    c$selected_edge_count <- retained_degree
    c$selected_removed_edge_count <- degree - retained_degree
    c$selected_retained_mass_fraction <- ratio(retained, matched)
    c$selected_retained_assigned_fraction <- ratio(retained, own)
    r <- s <- numeric(n)
    for (i in seq_len(n)) {
        resolved <- !is.na(labels[i]) & !is.na(opposite)
        same <- resolved & !is.na(labels[i] == opposite) & labels[i] == opposite
        r[i] <- sum(selected[i, resolved])
        s[i] <- sum(selected[i, same])
    }
    total <- ratio(s, own)
    conditional <- ratio(s, r)
    fields <- list(resolved_mass = r, same_state_mass = s,
        different_state_mass = r - s, unresolved_label_mass = retained - r,
        conditional_stability = conditional,
        same_state_fraction_of_matched_mass = ratio(s, retained),
        total_mass_stability = total,
        endpoint_unavailable_cells = as.integer(is.na(labels)),
        total_estimable_cells = as.integer(own > 0),
        total_nonestimable_cells = as.integer(own == 0),
        conditional_estimable_cells = as.integer(r > 0),
        conditional_nonestimable_cells = as.integer(r == 0),
        mean_cell_total_mass_stability = total,
        mean_cell_conditional_stability = conditional)
    for (field in names(fields)) c[[paste0("selected_", field)]] <-
        if (has_states) fields[[field]] else rep(NA_real_, n)
    if (has_states) c$state <- labels
    c
}

test_that("support public tables preserve exact code-built schema, types and values", {
    w <- matrix(0L, 6L, 7L)
    w[1:2, 1:2] <- 1L
    w[3L, 3:4] <- c(6L, 4L)
    w[4:5, 5L] <- c(3L, 2L)
    cases <- list(
        list(w = w, groups = list(c(1L, 1L, 2L, 3L, 3L, 4L),
                                  c(1L, 1L, 2L, 2L, 3L, 5L, 6L)),
             classes = c("many_to_many", "split", "merge", rep("unmatched", 3L)),
             extra = list(i = c(3L, 3L, NA_integer_), j = c(NA_integer_, NA_integer_, 5L),
                 s1 = c("assigned", "assigned", "unassigned"),
                 s2 = c("unassigned", "unassigned", "assigned"))),
        list(w = matrix(0L, 2L, 1L), groups = list(1:2, 3L),
             classes = rep("unmatched", 3L), extra = NULL),
        list(w = matrix(integer(), 0L, 0L), groups = list(integer(), integer()),
             classes = character(), extra = NULL),
        list(w = matrix(integer(), 0L, 2L), groups = list(integer(), 1:2),
             classes = rep("unmatched", 2L), extra = NULL))
    thresholds <- c(0, .01, .025, .05, .1, .2)
    pooled <- c(original_edges = "original_edges", original_matched_mass = "original_matched_mass",
        selected_edge_count = "retained_edges", selected_removed_edge_count = "removed_edges",
        selected_retained_mass = "retained_mass", selected_removed_mass = "excluded_matched_mass",
        selected_retained_mass_fraction = "retained_mass_fraction",
        selected_resolved_mass = "resolved_mass", selected_same_state_mass = "same_state_mass",
        selected_different_state_mass = "different_state_mass",
        selected_unresolved_label_mass = "unresolved_label_mass",
        selected_conditional_stability = "conditional_stability",
        selected_same_state_fraction_of_matched_mass = "same_state_fraction_of_matched_mass")
    for (case in cases) {
        w <- case$w
        cmp <- wp04_compare(w, extra = case$extra)
        own <- list(as.numeric(rowSums(w)), as.numeric(colSums(w)))
        if (!is.null(case$extra)) {
            own[[1]][3L] <- own[[1]][3L] + 2
            own[[2]][5L] <- own[[2]][5L] + 1
        }
        ids <- list(sprintf("x%d", seq_len(nrow(w))), sprintf("a%d", seq_len(ncol(w))))
        labs <- list(rep(c("s", "t", NA_character_), length.out = nrow(w)),
                     rep(c("s", NA_character_, "t"), length.out = ncol(w)))
        for (has in c(FALSE, TRUE)) {
            states <- if (has) stats::setNames(lapply(1:2, function(s)
                S4Vectors::DataFrame(cell_id = ids[[s]], state = labs[[s]])),
                c("left", "right")) else NULL
            before <- serialize(cmp, NULL)
            result <- supportSensitivity(cmp, states = states)
            expect_s4_class(result, "DataFrame")
            expect_identical(nrow(result), 12L)
            expect_identical(result$direction, rep(c("left", "right"), each = 6L))
            expect_identical(result$threshold, rep(thresholds, 2L))
            for (row in seq_len(nrow(result))) {
                side <- (row - 1L) %/% 6L + 1L
                threshold <- thresholds[(row - 1L) %% 6L + 1L]
                selected <- w
                for (i in seq_len(nrow(w))) for (j in seq_len(ncol(w))) {
                    if (w[i, j] < threshold * own[[side]][if (side == 1L) i else j])
                        selected[i, j] <- 0L
                }
                cells <- lapply(1:2, function(s) corrective_support_cells(
                    if (s == 1L) w else t(w), if (s == 1L) selected else t(selected),
                    own[[s]], labs[[s]], labs[[3L - s]], case$groups[[s]], ids[[s]], has))
                summaries <- lapply(cells, corrective_support_summary, has_states = has)
                scalar <- summaries[[side]]
                expect_identical(colnames(result), c("direction", "threshold", "states_supplied",
                                                     names(scalar), "cells", "components"))
                expect_identical(result$states_supplied[row], has)
                for (field in names(scalar)) expect_identical(result[[field]][row], scalar[[field]])
                actual <- result$cells[[row]]
                expect_s4_class(actual, "SimpleList")
                expect_identical(names(actual), c("left", "right"))
                for (s in 1:2) {
                    expect_identical(colnames(actual[[s]]), names(cells[[s]]))
                    for (field in names(cells[[s]])) expect_identical(actual[[s]][[field]], cells[[s]][[field]])
                }
                expected <- list(component_id = seq_along(case$classes),
                    n_segmentation_1 = tabulate(case$groups[[1]], length(case$classes)),
                    n_segmentation_2 = tabulate(case$groups[[2]], length(case$classes)),
                    component_class = case$classes)
                for (s in 1:2) {
                    groups <- lapply(seq_along(case$classes), function(k)
                        corrective_support_summary(cells[[s]][case$groups[[s]] == k, , drop = FALSE], has))
                    for (field in names(summaries[[s]])) {
                        expected[[paste0("selected_", field, "_", s)]] <-
                            vapply(groups, `[[`, summaries[[s]][[field]], field)
                    }
                    for (field in c("assigned_mass", "unmatched_mass"))
                        expected[[paste0(field, "_", s)]] <- vapply(groups, `[[`, 0, field)
                    if (s == 1L) for (key in names(pooled))
                        expected[[key]] <- vapply(groups, `[[`, scalar[[pooled[[key]]]], pooled[[key]])
                }
                components <- result$components[[row]]
                expect_s4_class(components, "DataFrame")
                expect_identical(nrow(components), length(case$classes))
                expect_identical(colnames(components), names(expected))
                for (field in names(expected)) expect_identical(components[[field]], unname(expected[[field]]))
            }
            expect_identical(serialize(cmp, NULL), before)
        }
    }
})
