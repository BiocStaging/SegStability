# Assisted-by: OpenAI Codex
test_that("absent geometry and unavailable sf remain explicitly unchecked", {
    cells <- S4Vectors::DataFrame(cell_id = "a")
    expect_identical(validation_row(validateGeometry(NULL, cells),
                                   "geometry.presence")$status, "NOT_APPLICABLE")
    testthat::local_mocked_bindings(.geometryAvailable = function() FALSE)
    r <- validateGeometry(structure(list(cell_id = "a"), class = "sf"), cells,
                          strict = TRUE)
    expect_identical(validation_row(r, "geometry.dependency")$status, "UNKNOWN")
    expect_identical(validation_row(r, "geometry.dependency")$severity, "warning")
})

validation_polygons <- function() {
    ring <- rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0))
    sf::st_sf(cell_id = c("a", "b"),
              geometry = sf::st_sfc(sf::st_polygon(list(ring)),
                                    sf::st_polygon(list(ring + 3))))
}

test_that("geometry warning failures never throw in strict mode or repair input", {
    skip_if_not_installed("sf")
    g <- validation_polygons()
    cells <- S4Vectors::DataFrame(cell_id = c("a", "b", "zero"))
    expect_false(validation_defect(validateGeometry(g, cells, strict = TRUE)))
    expect_identical(validation_row(validateGeometry(g, cells),
                                   "geometry.cell[3].presence")$status, "UNKNOWN")
    shapes <- list(sf::st_polygon(),
        sf::st_polygon(list(rbind(c(0, 0), c(2, 2), c(0, 2), c(2, 0), c(0, 0)))),
        sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(2, 0), c(0, 0)))))
    for (shape in shapes) {
        bad <- g
        sf::st_geometry(bad)[[1L]] <- shape
        before <- serialize(bad, NULL)
        r <- validateGeometry(bad, cells, strict = TRUE)
        expect_true(any(r$status == "FAIL" & r$severity == "warning"))
        expect_false(validation_defect(r))
        expect_identical(serialize(bad, NULL), before)
    }
    bad <- g
    bad$cell_id[1L] <- "unknown"
    expect_error(validateGeometry(bad, cells, strict = TRUE), "Validation failed")
    bad <- g
    bad$cell_id[2L] <- "a"
    expect_true(validation_defect(validateGeometry(bad, cells)))
    expect_identical(serialize(g, NULL), serialize(validation_polygons(), NULL))
})

test_that("holes multipart overlap zero populations and nonfinite geometry are explicit", {
    skip_if_not_installed("sf")
    g <- validation_polygons()
    cells <- S4Vectors::DataFrame(cell_id = c("a", "b"))
    outer <- unclass(sf::st_geometry(g)[[1L]])[[1L]]
    inner <- rbind(c(.5, .5), c(.5, 1.5), c(1.5, 1.5), c(1.5, .5), c(.5, .5))
    sf::st_geometry(g)[[1L]] <- sf::st_polygon(list(outer, inner))
    sf::st_geometry(g)[[2L]] <- sf::st_multipolygon(list(list(outer + 3),
                                                      list(outer + 6)))
    expect_false(validation_defect(validateGeometry(g[c(2L, 1L), ], cells)))
    bad <- validation_polygons()
    coords <- outer
    coords[2L, 1L] <- Inf
    sf::st_geometry(bad)[[1L]] <- sf::st_polygon(list(coords))
    before <- serialize(bad, NULL)
    r <- validateGeometry(bad, cells, strict = TRUE)
    expect_identical(validation_row(r, "geometry.cell[1].finite")$status, "FAIL")
    expect_identical(validation_row(r, "geometry.cell[1].finite")$severity, "warning")
    expect_identical(serialize(bad, NULL), before)
    expect_false(validation_defect(validateGeometry(g[integer(), ], cells)))
    empty_cells <- S4Vectors::DataFrame(cell_id = character())
    expect_false(validation_defect(validateGeometry(g[integer(), ], empty_cells)))
    expect_true(validation_defect(validateGeometry(g, empty_cells)))
    expect_true(validation_defect(validateGeometry(structure(list(), class = "sf"), cells)))
    points <- sf::st_sf(cell_id = c("a", "b"),
                       geometry = sf::st_sfc(sf::st_point(c(1, 1)),
                                             sf::st_point(c(2, 2))))
    expect_true(validation_defect(validateGeometry(points, cells)))
})

test_that("geographic local-unit conflicts fail without metadata repair", {
    skip_if_not_installed("sf")
    g <- validation_polygons()
    sf::st_crs(g) <- 4326
    before <- serialize(g, NULL)
    ci <- S4Vectors::SimpleList(units = "micrometres", crs = "EPSG:4326")
    expect_error(validateGeometry(g, S4Vectors::DataFrame(cell_id = c("a", "b")),
                                 ci, strict = TRUE), "Validation failed")
    expect_identical(serialize(g, NULL), before)
})

test_that("failed geometry prerequisites mark eligibility unchecked", {
    skip_if_not_installed("sf")
    g <- validation_polygons()
    g$cell_id[2L] <- "a"
    r <- validateGeometry(g, S4Vectors::DataFrame(cell_id = "a"))
    expect_identical(validation_row(r, "geometry.eligibility")$status, "UNKNOWN")
})


test_that("unverified backend validity cannot open a dependent area gate", {
    skip_if_not_installed("sf")
    g <- validation_polygons()
    local_mocked_bindings(st_is_valid = function(...) {
        warning("backend validity evidence unavailable")
        "Valid Geometry"
    }, .package = "sf")
    r <- validateGeometry(g, S4Vectors::DataFrame(cell_id = c("a", "b")),
                          strict = TRUE)
    expect_identical(validation_row(r, "geometry.cell[1].valid")$status, "UNKNOWN")
    expect_identical(validation_row(r, "geometry.cell[1].area")$status, "UNKNOWN")
    expect_false(.validationOK(validation_row(r, "geometry.cell[1].valid")))
})

test_that("reference ownership does not infer external geometry specimen", {
    skip_if_not_installed("sf")
    g <- validation_polygons()
    x <- wp03_set()
    r <- validateGeometry(g, S4Vectors::DataFrame(cell_id = c("a", "b")),
                          reference = x)
    expect_identical(validation_row(r, "geometry.specimen_id")$status, "UNKNOWN")
    r <- validateGeometry(g, S4Vectors::DataFrame(cell_id = c("a", "b")),
                          coordinate_info = S4Vectors::SimpleList(specimen_id = "other"),
                          reference = x)
    expect_error(.validationFinish(r, TRUE), "Validation failed")
})
