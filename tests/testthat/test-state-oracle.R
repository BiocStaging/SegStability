# Assisted-by: OpenAI Codex
test_that("historical portability missing endpoint mass stays unresolved", {
    # 144 portability: two 60/40 edges, one unavailable native endpoint.
    r <- wp06_result(matrix(c(60L, 40L), 2, 1), c("s", NA_character_), "s")
    expect_identical(r$global$matched_mass, 100)
    expect_identical(r$global$resolved_mass, 60)
    expect_identical(r$global$same_state_mass, 60)
    expect_identical(r$global$unresolved_label_mass, 40)
    expect_identical(r$global$different_state_mass, 0)
    expect_identical(r$global$conditional_stability, 1)
    expect_identical(r$global$total_mass_stability_1, .6)
    wp06_conservation(r)
})

test_that("historical estimand leakage partitions unmatched and unresolved", {
    # 121 estimand receipt: N=100, M=80, R=70, S=50.
    extra <- list(i = rep(1L, 20), j = rep(NA_integer_, 20),
                  s1 = rep("assigned", 20), s2 = rep("unassigned", 20))
    r <- wp06_result(matrix(c(50L, 20L, 10L), 1), "s",
                     c("s", "other", NA_character_), extra = extra)
    g <- r$global
    expect_identical(g$assigned_mass_1, 100)
    expect_identical(g$matched_mass, 80)
    expect_identical(g$resolved_mass, 70)
    expect_identical(g$same_state_mass, 50)
    expect_identical(g$different_state_mass, 20)
    expect_identical(g$unresolved_label_mass, 10)
    expect_identical(g$unmatched_mass_1, 20)
    expect_identical(g$conditional_stability, 50 / 70)
    expect_identical(g$same_state_fraction_of_matched_mass, 50 / 80)
    expect_identical(g$total_mass_stability_1, .5)
    wp06_conservation(r)
})

test_that("historical identity and zero cells retain estimable populations", {
    r <- wp06_result(diag(c(100L, 50L)), c(1L, 2L), c(1L, 2L))
    expect_identical(r$global$conditional_stability, 1)
    expect_identical(r$directional$mean_cell_total_mass_stability, c(1, 1))
    extra <- list(i = rep(2L, 20), j = rep(NA_integer_, 20),
                  s1 = rep("assigned", 20), s2 = rep("unassigned", 20))
    r <- wp06_result(matrix(0L, 2, 1), rep(NA_character_, 2), NA_character_,
                     extra = extra)
    d <- r$directional[1, ]
    expect_identical(d$n_cells, 2L)
    expect_identical(d$zero_molecule_cells, 1L)
    expect_identical(d$zero_correspondence_cells, 2L)
    expect_identical(d$endpoint_unavailable_cells, 2L)
    expect_identical(d$conditional_nonestimable_cells, 2L)
    expect_identical(d$total_mass_stability, 0)
    expect_identical(d$mean_cell_total_mass_stability, 0)
    expect_true(is.na(d$conditional_stability))
    wp06_conservation(r)
})
