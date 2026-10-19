# Assisted-by: OpenAI Codex
test_that("missing and ineligible geometry retain all cells and molecule edges", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(diag(2L, 3))
    g <- wp08_geometry(cmp)
    g$left <- g$left[1:2, ]
    g$right$geometry[[2]] <- sf::st_polygon()
    out <- wp09_run(cmp, g)
    expect_identical(out$cells$left$geometry_status,
                     c("eligible", "eligible", "unavailable"))
    expect_identical(out$cells$right$geometry_status,
                     c("eligible", "ineligible", "eligible"))
    expect_identical(length(out$edges), length(cmp@edges))
    expect_identical(nrow(out$cells$left), 3L)
    expect_true(all(is.na(out$pairs$iou[out$pairs$geometry_status != "eligible"])))
    expect_identical(out$summary$coverage_complete, FALSE)
    g <- list(left = NULL, right = NULL)
    out <- wp09_run(cmp, g)
    expect_identical(nrow(out$pairs), 3L)
    expect_identical(nrow(out$matching), 0L)
    expect_true(all(is.na(out$pairs$iou)))
})

test_that("geometry identity and coordinate failures block direct analysis", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(diag(2L, 2))
    g <- wp08_geometry(cmp)
    expect_error(compareGeometry(cmp, g), "coordinate|declaration")
    bad <- g
    bad$left$cell_id[1] <- "unknown"
    expect_error(wp09_run(cmp, bad), "Validation failed")
    bad <- g
    bad$left$cell_id[2] <- bad$left$cell_id[1]
    expect_error(wp09_run(cmp, bad), "Validation failed")
    bad <- g
    bad$left$cell_id <- 1:2
    expect_error(wp09_run(cmp, bad), "Validation failed")
    expect_error(wp09_run(cmp, unname(g)), "named list")
    info <- wp08_info()
    info$specimen_id <- "different"
    expect_error(compareGeometry(cmp, g, info), "conflict|prerequisite")
    bad <- g
    sf::st_crs(bad$left) <- 4326
    expect_error(wp09_run(cmp, bad), "CRS|Validation failed")
    expect_error(compareGeometry(NULL, g), "prerequisite")
})

test_that("holes and every multipolygon part contribute without repair", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(1L, 1, 1))
    outer <- unclass(wp08_box(0, 4, 0, 4))[[1L]]
    inner <- unclass(wp08_box(1, 3, 1, 3))[[1L]]
    hole <- sf::st_polygon(list(outer, inner))
    g <- wp09_geometry(cmp, list(hole), list(wp08_box(1, 3, 1, 3)))
    out <- wp09_run(cmp, g)
    expect_equal(out$pairs$area_1, 12)
    expect_equal(out$pairs$intersection_area, 0)
    mp <- sf::st_multipolygon(list(unclass(wp08_box(0, 1, 0, 1)),
                                  unclass(wp08_box(3, 4, 0, 1))))
    g <- wp09_geometry(cmp, list(mp), list(wp08_box(0, 4, 0, 1)))
    expect_equal(wp09_run(cmp, g)$pairs$intersection_area, 2)
    bow <- sf::st_polygon(list(rbind(c(0, 0), c(2, 2), c(0, 2), c(2, 0), c(0, 0))))
    g <- wp09_geometry(cmp, list(bow), list(wp08_box()))
    before <- serialize(g, NULL)
    expect_identical(wp09_run(cmp, g)$cells$left$geometry_status, "ineligible")
    expect_identical(serialize(g, NULL), before)
})

test_that("sf absence fails spatial entrypoint and leaves core operational", {
    cmp <- wp04_compare(matrix(1L, 1, 1))
    testthat::local_mocked_bindings(.geometryAvailable = function() FALSE,
                                  .package = "SegStability")
    expect_error(compareGeometry(cmp, list(left = NULL, right = NULL), wp08_info()),
                 "optional sf")
    expect_identical(length(correspondenceEdges(cmp)), 1L)
})

test_that("empty censuses and zero-edge populations are valid geometry inputs", {
    skip_if_not_installed("sf")
    for (mat in list(matrix(integer(), 0, 0), matrix(integer(), 2, 0),
                     matrix(integer(), 0, 2), matrix(0L, 2, 2))) {
        cmp <- wp04_compare(mat)
        g <- wp08_geometry(cmp)
        out <- wp09_run(cmp, g)
        expect_identical(out$components, cmp@components)
        expect_identical(length(out$edges), 0L)
        expect_identical(nrow(out$cells$left), nrow(mat))
        expect_identical(nrow(out$cells$right), ncol(mat))
        expect_identical(nrow(out$matching), min(nrow(mat), ncol(mat)))
    }
})

test_that("planar declarations and geometry CRS are checked without inference", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(1L, 1, 1))
    g <- wp08_geometry(cmp)
    info <- wp08_info()
    for (field in c("units", "origin", "axis_directions", "transforms",
                    "registration", "specimen_id")) {
        missing <- info
        missing[[field]] <- NULL
        expect_error(compareGeometry(cmp, g, missing), "declaration|prerequisite")
    }
    cmp@coordinate_info <- info
    changed <- info
    changed$origin <- c(1, 0)
    expect_error(compareGeometry(cmp, g, changed), "conflict")
    sf::st_crs(g$left) <- 3857
    sf::st_crs(g$right) <- 3857
    info$crs <- "EPSG:3857"
    expect_equal(compareGeometry(cmp, g, info)$pairs$iou, 1)
    sf::st_crs(g$right) <- NA
    expect_error(compareGeometry(cmp, g, info), "Validation failed|CRS")
    sf::st_crs(g$left) <- NA
    sf::st_crs(g$left) <- 4326
    sf::st_crs(g$right) <- 4326
    info$units <- "degrees"
    info$crs <- "EPSG:4326"
    cmp@coordinate_info <- S4Vectors::SimpleList()
    expect_error(compareGeometry(cmp, g, info), "Geographic|planar")
})

test_that("nonfinite nonpolygonal and unchecked geometry cannot confer metrics", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(1L, 1, 1))
    g <- wp08_geometry(cmp)
    coords <- unclass(wp08_box())[[1L]]
    coords[2, 1] <- Inf
    sf::st_geometry(g$left)[[1]] <- sf::st_polygon(list(coords))
    out <- wp09_run(cmp, g)
    expect_identical(out$cells$left$geometry_status, "ineligible")
    expect_true(is.na(out$pairs$iou))
    sf::st_geometry(g$left)[[1]] <- sf::st_point(c(1, 1))
    expect_error(wp09_run(cmp, g), "Validation failed")
    g <- wp08_geometry(cmp)
    testthat::local_mocked_bindings(st_is_valid = function(...) {
        warning("Validity unavailable")
        "Valid Geometry"
    }, .package = "sf")
    out <- wp09_run(cmp, g)
    expect_identical(out$cells$left$geometry_status, "unchecked")
    expect_identical(out$summary$matching_status, "unchecked")
    expect_true(is.na(out$pairs$iou))
})
