test_that("all five categorical endpoint cases are distinct", {
    labels <- list(c("s", "s"), c("s", "t"), c("s", NA_character_),
                   c(NA_character_, "s"), rep(NA_character_, 2))
    expected <- list(c(7, 7, 0, 0), c(7, 0, 7, 0), c(0, 0, 0, 7),
                     c(0, 0, 0, 7), c(0, 0, 0, 7))
    for (i in seq_along(labels)) {
        r <- wp06_result(matrix(7L), labels[[i]][1], labels[[i]][2])
        fields <- c("resolved_mass", "same_state_mass", "different_state_mass",
                    "unresolved_label_mass")
        expect_identical(unlist(as.list(r$global[, fields]), use.names = FALSE),
                         expected[[i]])
        e <- S4Vectors::mcols(r$edges)
        expect_identical(e$resolved, i <= 2L)
        expect_identical(e$same_state, if (i <= 2L) i == 1L else NA)
        wp06_conservation(r)
    }
})

test_that("splits merges and many-to-many retain mixed numeric evidence", {
    for (transpose in c(FALSE, TRUE)) {
        m <- matrix(c(6L, 4L), 1)
        if (transpose) m <- t(m)
        for (other in list(c("s", "s"), c("s", "t"),
                           c("s", NA_character_), rep(NA_character_, 2))) {
            r <- if (transpose) wp06_result(m, other, "s") else
                wp06_result(m, "s", other)
            expect_identical(r$global$same_state_mass,
                             if (all(is.na(other))) 0 else
                                 if (identical(other, c("s", "s"))) 10 else 6)
            expect_identical(r$components$component_class,
                             if (transpose) "merge" else "split")
            expect_identical(r$components$assigned_mass_1,
                             r$global$assigned_mass_1)
            wp06_conservation(r)
        }
    }
    r <- wp06_result(matrix(c(5L, 2L, 3L, 4L), 2), c("s", "t"), c("s", "t"))
    expect_identical(r$global$same_state_mass, 9)
    expect_identical(r$global$different_state_mass, 5)
    expect_identical(r$components$component_class, "many_to_many")
    wp06_conservation(r)
})

test_that("zero denominators and disconnected census cells stay visible", {
    r <- wp06_result(matrix(c(3L, 0L, 0L, 0L), 2), c("s", "unused"),
                     c("s", NA_character_))
    expect_identical(r$cells$left$matched_mass, c(3, 0))
    expect_identical(r$cells$left$total_mass_stability, c(1, NA_real_))
    expect_identical(r$components$matched_mass, c(3, 0, 0))
    expect_identical(r$directional$total_estimable_cells, c(1L, 1L))
    expect_identical(r$directional$mean_cell_conditional_stability, c(1, 1))
    wp06_conservation(r)
    r <- wp06_result(matrix(integer(), 0, 0), character(), character())
    expect_identical(nrow(r$cells$left), 0L)
    expect_identical(nrow(r$components), 0L)
    expect_identical(length(r$edges), 0L)
    expect_true(is.na(r$global$conditional_stability))
    expect_true(all(is.na(r$directional$mean_cell_total_mass_stability)))
    wp06_conservation(r)
})

test_that("supplied state types survive internal text comparison", {
    cmp <- wp04_compare(diag(c(3L, 2L, 1L)))
    l <- factor(c("1", "01", NA), levels = c("unused", "01", "1"), ordered = TRUE)
    rr <- c(1L, 1L, NA_integer_)
    states <- wp06_states(cmp, l, rr)
    r <- stateStability(cmp, states)
    expect_identical(r$cells$left$state, l)
    expect_identical(r$cells$right$state, rr)
    expect_identical(S4Vectors::mcols(r$edges)$state_1, l)
    expect_identical(S4Vectors::mcols(r$edges)$state_2, rr)
    expect_identical(r$global$same_state_mass, 3)
    expect_identical(r$global$different_state_mass, 2)
    expect_identical(r$global$unresolved_label_mass, 1)
    text <- c("1", "01", NA_character_)
    r <- stateStability(cmp, wp06_states(cmp, text, rr))
    expect_identical(r$cells$left$state, text)
    expect_identical(S4Vectors::mcols(r$edges)$state_1, text)
    # Different factor codes for the same labels must not alter equality.
    r <- stateStability(cmp, wp06_states(cmp, l,
        factor(c("1", "01", NA), levels = c("1", "01"))))
    expect_identical(r$global$same_state_mass, 5)
})

test_that("ambiguity union counts each molecule once", {
    statuses <- c("assigned", "unassigned", "ambiguous")
    s1 <- rep(statuses, each = 3L)
    s2 <- rep(statuses, 3L)
    extra <- list(i = ifelse(s1 == "assigned", 1L, NA_integer_),
                  j = ifelse(s2 == "assigned", 1L, NA_integer_), s1 = s1, s2 = s2)
    r <- wp06_result(matrix(2L), "s", "s", extra = extra)
    expect_identical(r$global$ambiguous_mass_1, 3)
    expect_identical(r$global$ambiguous_mass_2, 3)
    expect_identical(r$global$ambiguous_mass_union, 5)
    expect_identical(r$global$unmatched_mass_1, 2)
    expect_identical(r$global$assigned_mass_1, 5)
    wp06_conservation(r)
})

test_that("unavailable source evidence is UNKNOWN and nonblocking", {
    cmp <- wp04_compare()
    states <- wp06_states(cmp, c("s", "t"), c("s", "t"))
    r <- stateStability(cmp, states)
    report <- r$validation
    unavailable <- c("source.fingerprint_recomputation", "source.census_completeness",
                     "provenance.authenticity")
    expect_true(all(unavailable %in% report$check))
    expect_true(all(report$status[match(unavailable, report$check)] == "UNKNOWN"))
    expect_true(all(report$status[startsWith(report$check, "state.")] == "PASS"))
    expect_identical(r$global$same_state_mass, 5)
    bad <- cmp
    bad@cell_summary_1$assigned_mass[1] <- 1
    expect_error(stateStability(bad, states), "prerequisite|Validation")
    # WP05 warning capture yields UNKNOWN, which must close a required gate.
    dominance <- SegStability:::.directionalDominance
    testthat::local_mocked_bindings(.directionalDominance = function(...) {
        warning("unavailable prerequisite")
        dominance(...)
    }, .package = "SegStability")
    expect_error(stateStability(cmp, states), "UNKNOWN|prerequisite")
})

test_that("required state alignment warnings also block analysis", {
    cmp <- wp04_compare()
    states <- wp06_states(cmp, c("s", "t"), c("s", "t"))
    align <- SegStability:::.stateAlign
    testthat::local_mocked_bindings(.stateAlign = function(...) {
        warning("alignment evidence unavailable")
        align(...)
    }, .package = "SegStability")
    expect_error(stateStability(cmp, states), "UNKNOWN.*alignment|alignment.*UNKNOWN")
})

test_that("complete dominance and component evidence cannot be bypassed", {
    cmp <- wp04_compare()
    states <- wp06_states(cmp, c("s", "t"), c("s", "t"))
    bad <- cmp
    e <- S4Vectors::mcols(bad@edges)
    for (field in c("dominant_1", "dominant_2", "unique_dominant_1",
                    "unique_dominant_2", "reciprocal_dominant")) e[[field]][] <- FALSE
    S4Vectors::mcols(bad@edges) <- e
    bad@cell_summary_1$dominance_count[] <- 0L
    bad@cell_summary_2$dominance_count[] <- 0L
    expect_true(methods::validObject(bad, test = TRUE))
    expect_error(stateStability(bad, states), "dominance")
    bad <- cmp
    bad@cell_summary_1$component_id[] <- 1L
    bad@cell_summary_2$component_id[] <- 1L
    bad@components <- S4Vectors::DataFrame(component_id = 1L,
        n_segmentation_1 = 2L, n_segmentation_2 = 2L, component_class = "many_to_many")
    expect_true(methods::validObject(bad, test = TRUE))
    expect_error(stateStability(bad, states), "connected components")
})

test_that("equal-cell and mass-weighted estimands use explicit populations", {
    r <- wp06_result(diag(c(7L, 1L)), c("s", "t"), c("s", "s"))
    expect_identical(r$directional$total_mass_stability, c(7 / 8, 7 / 8))
    expect_identical(r$directional$conditional_stability, c(7 / 8, 7 / 8))
    expect_identical(r$directional$mean_cell_total_mass_stability, c(.5, .5))
    expect_identical(r$directional$mean_cell_conditional_stability, c(.5, .5))
    r <- wp06_result(diag(c(7L, 1L)), c("s", NA_character_), c("s", "s"))
    expect_identical(r$directional$mean_cell_total_mass_stability, c(.5, .5))
    expect_identical(r$directional$mean_cell_conditional_stability, c(1, 1))
    expect_identical(r$directional$conditional_estimable_cells, c(1L, 1L))
    expect_identical(r$directional$endpoint_unavailable_cells, c(1L, 0L))
})

test_that("cell and category permutations preserve scientific evidence", {
    m <- matrix(c(5L, 2L, 3L, 4L), 2)
    r <- wp06_result(m, c("s", "t"), c("s", "t"))
    renamed <- wp06_result(m, c("zebra", "alpha"), c("zebra", "alpha"),
                           ids1 = c("new1", "new2"), ids2 = c("other1", "other2"))
    permuted <- wp06_result(m[2:1, 2:1], c("t", "s"), c("t", "s"),
                            ids1 = c("x2", "x1"), ids2 = c("a2", "a1"))
    expect_identical(r$global, renamed$global)
    expect_identical(r$global, permuted$global)
    fields <- setdiff(colnames(r$cells$left), c("cell_id", "state"))
    expect_identical(r$cells$left[, fields], permuted$cells$left[2:1, fields])
    cmp <- wp04_compare()
    r <- stateStability(cmp, list(left = c(x1 = "", x2 = "NA"),
                                 right = c(a1 = "", a2 = NA_character_)))
    expect_identical(r$global$same_state_mass, 3)
    expect_identical(r$global$unresolved_label_mass, 2)
})

test_that("one empty census retains assigned unmatched cells on the other side", {
    extra <- list(i = c(1L, 1L), j = rep(NA_integer_, 2),
                  s1 = rep("assigned", 2), s2 = rep("unassigned", 2))
    r <- wp06_result(matrix(integer(), 2, 0), c("s", NA_character_), character(),
                     extra = extra)
    expect_identical(r$global$unmatched_mass_1, 2)
    expect_identical(r$global$assigned_mass_2, 0)
    expect_identical(r$components$assigned_mass_1, c(2, 0))
    expect_identical(r$components$assigned_mass_2, c(0, 0))
    expect_true(all(is.na(r$components$conditional_stability)))
    expect_identical(r$cells$left$total_mass_stability, c(0, NA_real_))
    wp06_conservation(r)
})

test_that("state mappings reject malformed identities shapes and values", {
    cmp <- wp04_compare()
    valid <- wp06_states(cmp, c("s", "t"), c("s", "t"))
    expect_error(stateStability(list(), valid), "SegmentationComparison|prerequisite")
    for (bad in list(NULL, unname(valid), valid[1], c(valid, list(extra = valid[[1]])),
                     stats::setNames(valid, c("left", "left")))) {
        expect_error(stateStability(cmp, bad))
    }
    badtables <- list(valid$left[1, ], valid$left[c(1, 1), ],
        S4Vectors::DataFrame(cell_id = c("x1", "unknown"), state = c("s", "t")),
        S4Vectors::DataFrame(cell_id = c("x1", NA_character_), state = c("s", "t")),
        S4Vectors::DataFrame(cell_id = c("x1", ""), state = c("s", "t")),
        S4Vectors::DataFrame(cell_id = 1:2, state = c("s", "t")),
        S4Vectors::DataFrame(cell_id = c("x1", "x2"), wrong = c("s", "t")))
    for (table in badtables) {
        bad <- valid
        bad$left <- table
        expect_error(stateStability(cmp, bad))
    }
    for (value in list(c(1, 2), c(TRUE, NA), list("s", "t"), as.Date(c("2020-01-01",
        "2020-01-02")), c(1 + 0i, 2 + 0i), as.raw(1:2), matrix(c("s", "t")))) {
        bad <- valid
        bad$left <- stats::setNames(value, c("x1", "x2"))
        expect_error(stateStability(cmp, bad))
    }
})

test_that("typed IDs are lossless and vector names never redefine their family", {
    families <- list(1:2, S4Vectors::as.LLint(c("9007199254740992", "9007199254740993")))
    if (requireNamespace("bit64", quietly = TRUE)) {
        families[[3]] <- bit64::as.integer64(c("9007199254740992", "9007199254740993"))
    }
    for (ids in families) {
        cmp <- wp04_compare(ids1 = ids, ids2 = rev(ids))
        states <- wp06_states(cmp, c("s", "t"), c("s", "t"))
        states$left <- states$left[2:1, ]
        r <- stateStability(cmp, states)
        expect_identical(r$cells$left$cell_id, ids)
        expect_identical(r$cells$left$state, c("s", "t"))
        states$left <- stats::setNames(c("s", "t"), as.character(ids))
        expect_error(stateStability(cmp, states), "character|family|typed")
    }
})

test_that("outputs follow comparison order and do not mutate or rematch", {
    cmp <- wp04_compare(matrix(c(3L, 1L, 2L, 4L), 2))
    before <- serialize(cmp, NULL)
    states <- list(right = c(a2 = "t", a1 = "s"), left = c(x2 = "t", x1 = "s"))
    r <- stateStability(cmp, states)
    expect_identical(names(r), c("global", "directional", "cells", "edges",
                                "components", "validation"))
    expect_s4_class(r, "SimpleList")
    expect_identical(names(r$cells), c("left", "right"))
    expect_identical(r$directional$segmentation, c("left", "right"))
    expect_identical(S4Vectors::from(r$edges), S4Vectors::from(correspondenceEdges(cmp)))
    expect_identical(S4Vectors::mcols(r$edges)[, colnames(S4Vectors::mcols(cmp@edges))],
                     S4Vectors::mcols(cmp@edges))
    expect_identical(serialize(cmp, NULL), before)
    expect_identical(stateStability(cmp, rev(states)), r)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(r, path)
    expect_identical(readRDS(path), r)
    # Equal run-local cell spellings do not create cross-run correspondence.
    same_ids <- wp04_compare(matrix(c(0L, 3L, 2L, 0L), 2),
                            ids1 = c("a", "b"), ids2 = c("a", "b"))
    r <- stateStability(same_ids, list(left = c(a = "s", b = "t"),
                                      right = c(a = "s", b = "t")))
    expect_identical(r$global$same_state_mass, 0)
    expect_identical(r$global$different_state_mass, 5)
    testthat::local_mocked_bindings(.sharedPairs = function(...) stop("rematch"),
        .aggregatePairs = function(...) stop("rematch"), .package = "SegStability")
    expect_no_error(stateStability(cmp, states))
})
