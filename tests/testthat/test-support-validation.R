# Assisted-by: OpenAI Codex
test_that("support rejects invalid thresholds and directions without coercion", {
    cmp <- wp04_compare()
    bad <- list(NULL, numeric(), TRUE, "0.1", NA_real_, NaN, Inf, -.1, 1.1,
                c(.1, .1), matrix(.1), factor("0.1"), structure(.1, class = "custom"))
    for (x in bad) expect_error(supportSensitivity(cmp, x), "threshold|prerequisite")
    for (x in list("", "LEFT", "unknown", NA_character_, 1L,
                   c("left", "right"), character(), factor("left"))) {
        expect_error(supportSensitivity(cmp, direction = x), "direction|prerequisite")
    }
    expect_identical(supportSensitivity(cmp, c(0L, 1L))$threshold, c(0, 1, 0, 1))
    expect_error(supportSensitivity(list()), "Comparison|comparison|prerequisite")
    expect_error(supportSensitivity(cmp, states = list()), "state|mapping|prerequisite")
})

test_that("scientific comparison corruption blocks support accounting", {
    cmp <- wp04_compare()
    for (field in c("weight", "fraction_1", "fraction_2", "dominant_1",
                    "unique_dominant_2", "reciprocal_dominant")) {
        bad <- cmp
        e <- S4Vectors::mcols(bad@edges)
        e[[field]][1] <- if (is.logical(e[[field]])) !e[[field]][1] else -1
        S4Vectors::mcols(bad@edges) <- e
        expect_error(supportSensitivity(bad), "prerequisite")
    }
    bad <- cmp
    bad@cell_summary_1$assigned_mass[1] <- 100
    expect_error(supportSensitivity(bad), "prerequisite")
    bad <- cmp
    e <- S4Vectors::mcols(bad@edges)
    e$duplicate <- e$weight
    colnames(e)[ncol(e)] <- "weight"
    S4Vectors::mcols(bad@edges) <- e
    expect_error(supportSensitivity(bad), "column|prerequisite")
})

test_that("only required UNKNOWN evidence blocks support", {
    cmp <- wp04_compare()
    real <- .directionalDominance
    testthat::local_mocked_bindings(.directionalDominance = function(...) {
        warning("Unavailable required dominance evidence")
        real(...)
    })
    expect_error(supportSensitivity(cmp), "UNKNOWN")
})

test_that("support uses no molecule rematching or stateStability call", {
    cmp <- wp04_compare()
    maps <- wp06_states(cmp, c("s", "t"), c("s", "t"))
    testthat::local_mocked_bindings(
        .sharedPairs = function(...) stop("Molecule rematching forbidden"),
        stateStability = function(...) stop("Nested stateStability forbidden"))
    expect_identical(supportSensitivity(cmp, states = maps)$retained_mass, rep(5, 12L))
})
