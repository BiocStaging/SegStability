# Assisted-by: OpenAI Codex
test_that("comparison scientific tables require unambiguous column names", {
    cmp <- wp04_compare()
    states <- wp06_states(cmp, c("s", "t"), c("s", "t"))
    bad <- cmp
    duplicate <- bad@cell_summary_1
    duplicate[["extra"]] <- duplicate$cell_id
    colnames(duplicate)[ncol(duplicate)] <- "cell_id"
    bad@cell_summary_1 <- duplicate
    expect_error(stateStability(bad, states), "column|prerequisite")
    bad <- cmp
    duplicate <- S4Vectors::mcols(bad@edges)
    duplicate[["extra"]] <- duplicate$weight
    colnames(duplicate)[ncol(duplicate)] <- "weight"
    S4Vectors::mcols(bad@edges) <- duplicate
    expect_error(stateStability(bad, states), "column|prerequisite")
})

test_that("complete tied graphs retain all edges and original dominance", {
    cmp <- wp04_compare(matrix(1L, 2, 2))
    states <- wp06_states(cmp, c("s", "t"), c("s", "t"))
    r <- stateStability(cmp, states)
    expect_identical(length(r$edges), 4L)
    expect_identical(r$cells$left$dominance_count, c(2L, 2L))
    expect_true(all(S4Vectors::mcols(r$edges)$dominant_1))
    expect_false(any(S4Vectors::mcols(r$edges)$reciprocal_dominant))
    expect_identical(r$global$same_state_mass, 2)
    expect_identical(r$cells$left$total_mass_stability, c(.5, .5))
    expect_identical(r$components$matched_mass, 4)
    wp06_conservation(r)
})

test_that("opaque spellings and invalid vector names are never repaired", {
    ids <- c("9223372036854775809", "9223372036854775810")
    cmp <- wp04_compare(ids1 = ids)
    states <- list(left = stats::setNames(c(" S ", "s"), ids),
                   right = c(a1 = "S", a2 = "s"))
    r <- stateStability(cmp, states)
    expect_identical(r$cells$left$cell_id, ids)
    expect_identical(r$global$same_state_mass, 2)
    expect_identical(r$global$different_state_mass, 3)
    for (keys in list(rep(ids[1], 2), c(ids[1], ""), c(ids[1], NA_character_),
                      c(ids[1], "unknown"), NULL)) {
        bad <- states
        bad$left <- stats::setNames(c("s", "t"), keys)
        expect_error(stateStability(cmp, bad))
    }
})

test_that("state validation and source limitations survive serialization", {
    cmp <- wp04_compare(matrix(c(3L, 2L), 1))
    states <- wp06_states(cmp, factor("s", levels = c("unused", "s")), c(1L, 2L))
    r <- stateStability(cmp, states)
    expect_identical(S4Vectors::mcols(r$edges)$state_1,
                     rep(states$left$state, 2))
    expect_identical(S4Vectors::mcols(r$edges)$state_2, states$right$state)
    expect_identical(names(formals(stateStability)), c("comparison", "states"))
    f <- tempfile(fileext = ".rds")
    on.exit(unlink(f))
    saveRDS(cmp, f)
    expect_identical(stateStability(readRDS(f), states), r)
    expect_false(any(r$validation$status[
        r$validation$check == "source.fingerprint_recomputation"] == "PASS"))
})

test_that("a factor NA label is unresolved even when its integer code exists", {
    cmp <- wp04_compare()
    value <- factor(c(NA, "s"), exclude = NULL)
    # R's is.na.factor sees codes; categorical availability must see labels.
    expect_false(anyNA(value))
    r <- stateStability(cmp, wp06_states(cmp, value, value))
    expect_identical(r$cells$left$state, value)
    expect_identical(S4Vectors::mcols(r$edges)$state_1, value)
    expect_identical(r$global$resolved_mass, 2)
    expect_identical(r$global$unresolved_label_mass, 3)
    expect_identical(r$global$different_state_mass, 0)
    expect_identical(r$directional$endpoint_unavailable_cells, c(1L, 1L))
    expect_identical(S4Vectors::mcols(r$edges)$same_state, c(NA, TRUE))
    wp06_conservation(r)
})
