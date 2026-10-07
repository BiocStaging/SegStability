test_that("whole polygon bboxes distinguish interior from touching crossing and outside", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(3, 1, 1))
    shapes <- list(wp08_box(), wp08_box(0, 2), wp08_box(-1, 2), wp08_box(11, 12))
    for (i in seq_along(shapes)) {
        g <- wp08_geometry(cmp)
        sf::st_geometry(g$left)[[1]] <- shapes[[i]]
        r <- anchorTopology(cmp, g, wp08_window(), wp08_info())
        expect_identical(r$cells$left$strict_interior, i == 1L)
        expect_identical(nrow(r$anchors), as.integer(i == 1L))
        expect_identical(r$cells$left$geometry_status, "eligible")
        if (i > 1L) expect_identical(r$cells$left$selection_stratum,
                                     "reciprocal_noninterior")
        wp08_unchanged(cmp, g, anchorTopology(cmp, g, wp08_window(), wp08_info()))
    }
    g <- wp08_geometry(cmp)
    sf::st_geometry(g$right)[[1]] <- wp08_box(1, 10)
    expect_identical(nrow(anchorTopology(cmp, g, wp08_window(), wp08_info())$anchors), 0L)
})

test_that("holes and all multipolygon parts determine unaltered geometry eligibility", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(3, 1, 1))
    g <- wp08_geometry(cmp)
    outer <- unclass(wp08_box(1, 5, 1, 5))[[1]]
    hole <- rbind(c(2, 2), c(2, 3), c(3, 3), c(3, 2), c(2, 2))
    sf::st_geometry(g$left)[[1]] <- sf::st_polygon(list(outer, hole))
    parts <- list(unclass(wp08_box()), unclass(wp08_box(6, 7)))
    sf::st_geometry(g$right)[[1]] <- sf::st_multipolygon(parts)
    expect_identical(anchorTopology(cmp, g, wp08_window(), wp08_info())$anchors$anchor_topology,
                     "one_to_one_like")
    parts[[2]] <- unclass(wp08_box(9, 11))
    sf::st_geometry(g$right)[[1]] <- sf::st_multipolygon(parts)
    expect_false(anchorTopology(cmp, g, wp08_window(), wp08_info())$cells$right$strict_interior)
    wp08_unchanged(cmp, g, anchorTopology(cmp, g, wp08_window(), wp08_info()))
})

test_that("invalid empty nonfinite zero area and missing geometry are explicit retained strata", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(3, 1, 1))
    invalid <- sf::st_polygon(list(rbind(c(1, 1), c(3, 3), c(1, 3), c(3, 1), c(1, 1))))
    coords <- unclass(wp08_box())[[1]]
    coords[2, 1] <- Inf
    shapes <- list(invalid, sf::st_polygon(), sf::st_polygon(list(coords)),
        sf::st_polygon(list(rbind(c(1, 1), c(2, 1), c(3, 1), c(1, 1)))))
    for (shape in shapes) {
        g <- wp08_geometry(cmp)
        sf::st_geometry(g$left)[[1]] <- shape
        before <- serialize(g, NULL)
        r <- anchorTopology(cmp, g, wp08_window(), wp08_info())
        expect_identical(r$cells$left$geometry_status, "ineligible")
        expect_false(r$cells$left$geometry_eligible)
        expect_identical(r$cells$left$selection_stratum, "reciprocal_geometry_ineligible")
        expect_identical(r$cells$right$selection_stratum, "reciprocal_geometry_ineligible")
        expect_identical(nrow(r$cells$left), 1L)
        expect_identical(length(r$edges), 1L)
        expect_identical(nrow(r$anchors), 0L)
        expect_identical(serialize(g, NULL), before)
    }
    g <- wp08_geometry(cmp)
    sf::st_geometry(g$left)[[1]] <- invalid
    anchorTopology(cmp, g, wp08_window(), wp08_info())
    expect_false(sf::st_is_valid(g$left))
    for (absent in list(NULL, g$left[integer(), ])) {
        g$left <- absent
        # List replacement must retain the explicit NULL segmentation entry.
        if (!"left" %in% names(g)) g["left"] <- list(NULL)
        r <- anchorTopology(cmp, g, wp08_window(), wp08_info())
        expect_identical(r$cells$left$geometry_status, "unavailable")
        expect_identical(r$cells$right$selection_stratum, "reciprocal_geometry_unavailable")
        expect_true(is.na(r$cells$left$strict_interior))
        expect_identical(nrow(r$cells$left), 1L)
    }
})

test_that("identity type and declared CRS conflicts fail rather than repair", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare()
    g <- wp08_geometry(cmp)
    for (ids in list(c("x1", "x1"), c("x1", "unknown"), 1:2, c("x1", NA))) {
        bad <- g
        bad$left$cell_id <- ids
        expect_error(anchorTopology(cmp, bad, wp08_window(), wp08_info()), "Validation failed")
    }
    bad <- g
    sf::st_geometry(bad$left)[[1]] <- sf::st_point(c(1, 1))
    expect_error(anchorTopology(cmp, bad, wp08_window(), wp08_info()), "Validation failed")
    expect_error(anchorTopology(cmp, unname(g), wp08_window(), wp08_info()), "geometry")
    expect_error(anchorTopology(cmp, g[1], wp08_window(), wp08_info()), "geometry")
    bad <- g
    sf::st_crs(bad$left) <- 4326
    before <- serialize(bad, NULL)
    info <- wp08_info()
    info$crs <- "EPSG:4326"
    expect_error(anchorTopology(cmp, bad, wp08_window(), info))
    expect_identical(serialize(bad, NULL), before)
    bad <- g
    sf::st_crs(bad$left) <- 3857
    expect_error(anchorTopology(cmp, bad, wp08_window(), wp08_info()), "prerequisite")
    info$crs <- "EPSG:3857"
    expect_error(anchorTopology(cmp, bad, wp08_window(), info))
    sf::st_crs(bad$right) <- 3857
    expect_identical(nrow(anchorTopology(cmp, bad, wp08_window(), info)$anchors), 2L)
})

test_that("unverified geometry backend gates cannot confer anchor eligibility", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare()
    g <- wp08_geometry(cmp)
    local_mocked_bindings(st_is_valid = function(...) {
        warning("backend unavailable")
        "Valid Geometry"
    }, .package = "sf")
    r <- anchorTopology(cmp, g, wp08_window(), wp08_info())
    expect_identical(nrow(r$anchors), 0L)
    expect_identical(r$cells$left$geometry_status, rep("unchecked", 2))
    expect_true(all(is.na(r$cells$left$geometry_eligible)))
})

test_that("all audited ID families align exactly without numeric conversion", {
    skip_if_not_installed("sf")
    families <- list(c("01", "1"), 1:2,
        S4Vectors::as.LLint(c("9007199254740992", "9007199254740993")))
    if (requireNamespace("bit64", quietly = TRUE)) {
        families[[4]] <- bit64::as.integer64(c("9007199254740992", "9007199254740993"))
    }
    for (ids in families) {
        cmp <- wp04_compare(ids1 = ids, ids2 = ids)
        g <- wp08_geometry(cmp)
        r <- anchorTopology(cmp, lapply(g, function(x) x[2:1, ]), wp08_window(), wp08_info())
        expect_identical(r$anchors$cell_id_1, ids)
        expect_identical(r$anchors$cell_id_2, ids)
        expect_identical(r$cells$left$cell_id, ids)
        # Same geometry, different run-local IDs; correspondence is molecule-based.
        other <- wp04_compare(ids1 = ids, ids2 = ids[2:1])
        r2 <- anchorTopology(other, wp08_geometry(other), wp08_window(), wp08_info())
        expect_identical(r2$anchors$cell_id_2, ids[2:1])
        expect_identical(r2$anchors$anchor_topology, r$anchors$anchor_topology)
    }
})

test_that("topology never calls repair transformation clipping buffering or simplification", {
    skip_if_not_installed("sf")
    forbidden <- function(...) stop("Forbidden geometry alteration called")
    local_mocked_bindings(st_make_valid = forbidden, st_buffer = forbidden,
        st_intersection = forbidden, st_simplify = forbidden,
        st_transform = forbidden, .package = "sf")
    cmp <- wp04_compare(matrix(c(60, 40), 1))
    g <- wp08_geometry(cmp)
    sf::st_geometry(g$left)[[1]] <- wp08_box(-1, 11)
    expect_identical(nrow(anchorTopology(cmp, g, wp08_window(), wp08_info())$anchors), 0L)
    wp08_unchanged(cmp, g, anchorTopology(cmp, g, wp08_window(), wp08_info()))
})
