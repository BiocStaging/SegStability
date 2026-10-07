test_that("support signature, defaults and direction-major order are exact", {
    expect_identical(names(formals(supportSensitivity)),
                     c("comparison", "thresholds", "direction", "states"))
    cmp <- wp04_compare()
    out <- supportSensitivity(cmp)
    expect_s4_class(out, "DataFrame")
    expect_identical(out$direction, rep(c("left", "right"), each = 6L))
    expect_identical(out$threshold, rep(c(0, .01, .025, .05, .1, .2), 2L))
    expect_identical(out$retained_edges, rep(2L, 12L))
    expect_identical(out$retained_mass, rep(5, 12L))
    expect_true(all(is.na(out$resolved_mass)))
    expect_true(all(is.na(out$unresolved_label_mass)))
    expect_true(all(is.na(out$conditional_estimable_cells)))
    expect_identical(supportSensitivity(cmp, c(.5, 0, 1), "right")$threshold,
                     c(.5, 0, 1))
    report <- S4Vectors::metadata(out)$validation
    expect_true(all(report$status[startsWith(report$check, "support.")] == "PASS"))
    expect_identical(report$status[report$check == "source.fingerprint_recomputation"],
                     "UNKNOWN")
    expect_identical(report$status[report$check == "source.census_completeness"], "UNKNOWN")
})

test_that("directional fraction selection preserves independent denominators", {
    extra <- list(i = 1L, j = NA_integer_, s1 = "assigned", s2 = "unassigned")
    cmp <- wp04_compare(matrix(c(6L, 4L), 1L), extra = extra)
    before <- cmp
    out <- supportSensitivity(cmp, c(0, 4/11, 4/11 + 1e-12, 1),
                              states = wp06_states(cmp, "s", c("s", "t")))
    expect_identical(out$retained_edges, c(2L, 2L, 1L, 0L, rep(2L, 4L)))
    expect_identical(out$retained_mass, c(10, 10, 6, 0, rep(10, 4)))
    expect_identical(out$assigned_mass, c(rep(11, 4), rep(10, 4)))
    expect_identical(out$unmatched_mass, c(rep(1, 4), rep(0, 4)))
    expect_identical(out$retained_mass_fraction[3], .6)
    expect_identical(out$total_mass_stability[3], 6/11)
    expect_identical(out$conditional_stability[3], 1)
    expect_identical(out$cells[[3]]$left$assigned_mass, 11)
    expect_identical(out$cells[[3]]$left$matched_mass, 10)
    expect_identical(out$cells[[3]]$left$selected_retained_mass, 6)
    expect_identical(out$cells[[3]]$left$selected_retained_assigned_fraction, 6/11)
    expect_identical(cmp, before)
    for (i in seq_len(nrow(out))) wp07_conservation(out, i, cmp)
})

test_that("states cannot influence any structural support output", {
    cmp <- wp04_compare(matrix(c(5L, 2L, 1L, 3L, 0L, 0L), 2L),
                        ids1 = c("x", "z"), ids2 = c("a", "b", "empty"))
    masks <- list()
    real <- .supportCellTables
    testthat::local_mocked_bindings(.supportCellTables = function(
        comparison, weight, keep, labels) {
        masks[[length(masks) + 1L]] <<- keep
        real(comparison, weight, keep, labels)
    })
    a <- supportSensitivity(cmp, c(.7, 0, 1, .2))
    b <- supportSensitivity(cmp, c(.7, 0, 1, .2),
        states = wp06_states(cmp, c("s", NA_character_), c("s", "t", "u")))
    expect_identical(masks[seq_len(nrow(a))], masks[nrow(a) + seq_len(nrow(b))])
    state_fields <- c("states_supplied", wp07_state_fields())
    scalar <- setdiff(colnames(a), c(state_fields, "cells", "components"))
    # State mapping adds required validation rows, not structural evidence.
    structural_a <- a[, scalar]
    structural_b <- b[, scalar]
    S4Vectors::metadata(structural_a) <- list()
    S4Vectors::metadata(structural_b) <- list()
    expect_identical(structural_a, structural_b)
    for (i in seq_len(nrow(a))) {
        for (side in 1:2) {
            x <- a$cells[[i]][[side]]
            y <- b$cells[[i]][[side]]
            fields <- setdiff(colnames(x), c("state", paste0("selected_", wp07_state_fields())))
            expect_identical(x[, fields], y[, fields])
        }
        fields <- setdiff(colnames(a$components[[i]]), c(
            paste0("selected_", wp07_state_fields()),
            paste0("selected_", wp07_state_fields(), "_1"),
            paste0("selected_", wp07_state_fields(), "_2")))
        expect_identical(a$components[[i]][, fields], b$components[[i]][, fields])
    }
})

test_that("full censuses and original component identities survive removal", {
    cmp <- wp04_compare(matrix(c(5L, 1L, 0L, 4L, 0L, 0L), 2L))
    out <- supportSensitivity(cmp, 1)
    for (i in seq_len(nrow(out))) {
        for (side in 1:2) {
            original <- cellSummary(cmp, segmentationNames(cmp)[side])
            cells <- out$cells[[i]][[side]]
            expect_identical(cells[, colnames(original)], original)
        }
        original <- componentSummary(cmp)
        expect_identical(out$components[[i]][, colnames(original)], original)
        wp07_conservation(out, i, cmp)
    }
    tied <- wp04_compare(matrix(1L, 2L, 2L))
    before <- correspondenceEdges(tied)
    r <- supportSensitivity(tied, c(.5, .5 + 1e-12))
    expect_identical(r$retained_edges, c(4L, 0L, 4L, 0L))
    expect_identical(r$cells[[2]]$left$dominance_count, c(2L, 2L))
    expect_identical(r$components[[2]]$component_class, "many_to_many")
    expect_identical(correspondenceEdges(tied), before)
})

test_that("empty graphs and empty censuses remain explicitly nonestimable", {
    for (shape in list(c(0L, 0L), c(2L, 0L), c(0L, 2L), c(2L, 3L))) {
        cmp <- wp04_compare(matrix(0L, shape[1], shape[2]))
        out <- supportSensitivity(cmp, c(0, 1))
        expect_identical(out$retained_edges, rep(0L, 4L))
        expect_identical(out$retained_mass, rep(0, 4L))
        expect_true(all(is.na(out$retained_mass_fraction)))
        expect_true(all(is.na(out$retained_assigned_fraction)))
        for (i in seq_len(nrow(out))) wp07_conservation(out, i, cmp)
        maps <- wp06_states(cmp, rep("s", shape[1]), rep("s", shape[2]))
        r <- supportSensitivity(cmp, 1, states = maps)
        expect_identical(r$resolved_mass, c(0, 0))
        expect_true(all(is.na(r$conditional_stability)))
        expect_true(all(is.na(r$total_mass_stability)))
    }
})

test_that("NULL labels differ from explicitly unresolved labels", {
    cmp <- wp04_compare(matrix(c(6L, 4L), 1L))
    absent <- supportSensitivity(cmp, .5)
    unknown <- supportSensitivity(cmp, .5,
        states = wp06_states(cmp, "s", c(NA_character_, NA_character_)))
    expect_true(all(is.na(absent$unresolved_label_mass)))
    expect_identical(unknown$unresolved_label_mass, c(6, 10))
    expect_identical(unknown$resolved_mass, c(0, 0))
    expect_true(all(is.na(unknown$conditional_stability)))
    expect_identical(unknown$total_mass_stability, c(0, 0))
    values <- factor(c(NA, "s"), exclude = NULL)
    cmp <- wp04_compare()
    r <- supportSensitivity(cmp, 0, states = wp06_states(cmp, values, values))
    expect_identical(r$unresolved_label_mass, c(3, 3))
    expect_identical(r$cells[[1]]$left$state, values)
    expect_identical(r$cells[[1]]$left$selected_endpoint_unavailable_cells, c(1L, 0L))
})

test_that("unmatched assigned and ambiguous populations never acquire selected edges", {
    extra <- list(i = c(1L, NA_integer_, NA_integer_),
                  j = c(NA_integer_, 1L, NA_integer_),
                  s1 = c("assigned", "ambiguous", "ambiguous"),
                  s2 = c("ambiguous", "assigned", "ambiguous"))
    cmp <- wp04_compare(matrix(0L, 2L, 2L), extra = extra)
    status <- assignmentSummary(cmp)
    out <- supportSensitivity(cmp, c(0, 1),
        states = wp06_states(cmp, c("s", "t"), c("s", "t")))
    expect_identical(out$assigned_mass, rep(1, 4L))
    expect_identical(out$unmatched_mass, rep(1, 4L))
    expect_identical(out$total_mass_stability, rep(0, 4L))
    expect_identical(out$zero_correspondence_cells, rep(2L, 4L))
    expect_identical(out$zero_molecule_cells, rep(1L, 4L))
    expect_identical(assignmentSummary(cmp), status)
    for (i in seq_len(nrow(out))) wp07_conservation(out, i, cmp)
})
