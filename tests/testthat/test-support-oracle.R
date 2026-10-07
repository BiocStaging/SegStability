test_that("historical 60/40 support cost remains over original population", {
    cmp <- wp04_compare(matrix(c(60L, 40L), 1L))
    out <- supportSensitivity(cmp, .5, "left",
        wp06_states(cmp, "s", c("s", NA_character_)))
    expect_identical(out$retained_mass, 60)
    expect_identical(out$excluded_matched_mass, 40)
    expect_identical(out$assigned_mass, 100)
    expect_identical(out$total_mass_stability, .6)
    expect_identical(out$conditional_stability, 1)
    expect_identical(out$unresolved_label_mass, 0)
    expect_identical(out$components[[1]]$component_class, "split")
    wp07_conservation(out, 1, cmp)
})

test_that("historical structural matrices obey independently hand-computed retention", {
    matrices <- list(matrix(c(60L, 0L, 0L, 40L), 2L),
        matrix(c(60L, 40L), 1L), matrix(c(60L, 40L), 2L),
        matrix(c(50L, 20L, 30L, 0L), 2L), matrix(50L, 2L, 2L))
    left_mass <- c(100, 60, 100, 70, 0)
    right_mass <- c(100, 100, 60, 80, 0)
    for (i in seq_along(matrices)) {
        cmp <- wp04_compare(matrices[[i]])
        out <- supportSensitivity(cmp, .6)
        expect_identical(out$retained_mass, c(left_mass[i], right_mass[i]))
        for (j in 1:2) wp07_conservation(out, j, cmp)
    }
})

test_that("retention monotonicity uses threshold order, not component class assumptions", {
    cmp <- wp04_compare(matrix(c(5L, 1L, 2L, 4L), 2L))
    out <- supportSensitivity(cmp, c(0, .1, .2, .4, .6, 1))
    for (direction in c("left", "right")) {
        r <- out[out$direction == direction, ]
        expect_true(all(diff(r$retained_edges) <= 0))
        expect_true(all(diff(r$retained_mass) <= 0))
        expect_true(all(diff(r$removed_edges) >= 0))
        expect_true(all(diff(r$excluded_matched_mass) >= 0))
    }
    # Same mass can describe different original structure; no inferred ranking.
    a <- supportSensitivity(wp04_compare(matrix(c(5L, 5L), 1L)), .5, "left")
    b <- supportSensitivity(wp04_compare(diag(c(5L, 5L))), .5, "left")
    expect_identical(a$retained_mass, b$retained_mass)
    expect_identical(a$components[[1]]$component_class, "split")
    expect_identical(b$components[[1]]$component_class, c("one_to_one", "one_to_one"))
})

test_that("opaque ID families and serialization preserve deterministic results", {
    ids <- list(c("01", "1"), c(1L, 2L),
        S4Vectors::as.LLint(c("9007199254740992", "9007199254740993")),
        c("9223372036854775809", "9223372036854775810"))
    if (requireNamespace("bit64", quietly = TRUE)) {
        ids[[length(ids) + 1L]] <- bit64::as.integer64(c("9007199254740992",
                                                       "9007199254740993"))
    }
    for (x in ids) {
        cmp <- wp04_compare(ids1 = x)
        maps <- wp06_states(cmp, c(1L, 2L), factor(c("1", "2")))
        out <- supportSensitivity(cmp, c(1, 0), states = maps)
        expect_identical(out$cells[[1]]$left$cell_id, x)
        expect_identical(out$retained_mass, rep(5, 4L))
        expect_identical(out$same_state_mass, rep(5, 4L))
        f <- tempfile(fileext = ".rds")
        saveRDS(cmp, f)
        expect_identical(supportSensitivity(readRDS(f), c(1, 0), states = maps), out)
        unlink(f)
    }
})
