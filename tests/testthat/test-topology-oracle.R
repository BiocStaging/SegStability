test_that("historical reciprocal fan fixtures retain exact frozen class literals", {
    skip_if_not_installed("sf")
    # 136 validation.json and 144 portability oracle; no historical code copied.
    cases <- list(identity = matrix(c(60, 0, 0, 40), 2),
        split = matrix(c(60, 40), 1), merge = matrix(c(60, 40), 2),
        complex = matrix(c(50, 20, 30, 0), 2), tied = matrix(50, 2, 2))
    expected <- list(c("one_to_one_like", "one_to_one_like"), "split_like",
                     "merge_like", "complex", character())
    for (i in seq_along(cases)) {
        r <- wp08_run(cases[[i]])
        expect_identical(r$anchors$anchor_topology, expected[[i]])
        expect_type(S4Vectors::mcols(r$edges)$anchor_topology, "character")
        edge_data <- S4Vectors::mcols(r$edges)
        expect_identical(edge_data$anchor_topology[edge_data$interior_anchor], expected[[i]])
        for (side in 1:2) {
            expect_type(r$cells[[side]]$anchor_topology, "character")
            expect_identical(r$cells[[side]]$anchor_topology[r$cells[[side]]$interior_anchor],
                             expected[[i]])
        }
    }
    r <- wp08_run(cases$complex)
    expect_identical(r$anchors$fan_out, 2L)
    expect_identical(r$anchors$fan_in, 2L)
    expect_identical(S4Vectors::mcols(r$edges)$anchor_topology,
                     c("complex", "undefined", "undefined"))
    expect_identical(wp08_run(cases$tied)$cells[[1]]$selection_stratum,
                     rep("tied_dominance", 2))
    expect_identical(nrow(wp08_run(matrix(c(50, 50), 1))$anchors), 0L)
})

test_that("interior fan contributors need not be reciprocal and never redefine dominance", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(c(60, 40), 1))
    g <- wp08_geometry(cmp)
    r <- anchorTopology(cmp, g, wp08_window(), wp08_info())
    expect_identical(r$anchors$fan_out, 2L)
    sf::st_geometry(g$right)[[2]] <- wp08_box(0, 1)
    r <- anchorTopology(cmp, g, wp08_window(), wp08_info())
    expect_identical(r$anchors$anchor_topology, "one_to_one_like")
    g$right <- g$right[1, ]
    r2 <- anchorTopology(cmp, g, wp08_window(), wp08_info())
    expect_identical(r2$anchors$anchor_topology, "one_to_one_like")
    expect_identical(r2$cells$right$geometry_status, c("eligible", "unavailable"))
    for (field in colnames(S4Vectors::mcols(correspondenceEdges(cmp)))) {
        expect_identical(S4Vectors::mcols(r2$edges)[[field]],
                         S4Vectors::mcols(correspondenceEdges(cmp))[[field]])
    }
    # The dominant edge is boundary-excluded; the lesser edge cannot become an anchor.
    g <- wp08_geometry(cmp)
    sf::st_geometry(g$right)[[1]] <- wp08_box(0, 1)
    expect_identical(nrow(anchorTopology(cmp, g, wp08_window(), wp08_info())$anchors), 0L)
})

test_that("run-local relabeling and side reversal preserve directional fan definitions", {
    skip_if_not_installed("sf")
    m <- matrix(c(60, 40), 1)
    r <- wp08_run(m)
    relabeled <- wp08_run(m, ids1 = "other", ids2 = c("z", "a"))
    expect_identical(relabeled$anchors$anchor_topology, r$anchors$anchor_topology)
    expect_identical(relabeled$anchors$fan_out, r$anchors$fan_out)
    expect_identical(relabeled$anchors$fan_in, r$anchors$fan_in)
    cmp <- wp04_compare(t(m))
    reverse <- anchorTopology(cmp, wp08_geometry(cmp), wp08_window(), wp08_info())
    expect_identical(reverse$anchors$anchor_topology, "merge_like")
    expect_identical(reverse$anchors$fan_in, r$anchors$fan_out)
    expect_identical(reverse$anchors$fan_out, r$anchors$fan_in)
})
