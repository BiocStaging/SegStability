test_that("typed cell identities and row permutations are exact", {
    skip_if_not_installed("sf")
    ids <- list(1:2, c("001", "1"),
        S4Vectors::as.LLint(c("9007199254740992", "9007199254740993")))
    if (requireNamespace("bit64", quietly = TRUE)) {
        ids <- c(ids, list(bit64::as.integer64(c("9007199254740992",
                                               "9007199254740993"))))
    }
    for (id in ids) {
        cmp <- wp04_compare(diag(2L, 2), ids1 = id, ids2 = id)
        g <- wp08_geometry(cmp)
        out <- wp09_run(cmp, g)
        for (i in 1:2) {
            expect_identical(out$cells[[i]]$cell_id, id)
            expect_identical(out$matching[[paste0("cell_id_", i)]], id)
        }
        g <- lapply(g, function(x) x[2:1, ])
        shuffled <- wp09_run(cmp, g)
        expect_identical(out, shuffled)
        path <- tempfile(fileext = ".rds")
        on.exit(unlink(path), add = TRUE)
        saveRDS(out, path)
        expect_identical(readRDS(path), out)
    }
})

test_that("splits merges complex and tied graphs keep every molecule edge", {
    skip_if_not_installed("sf")
    for (mat in list(matrix(c(60L, 40L), 1), matrix(c(60L, 40L), 2),
                     matrix(c(50L, 20L, 30L, 0L), 2), matrix(50L, 2, 2))) {
        cmp <- wp04_compare(mat)
        g <- wp08_geometry(cmp)
        out <- wp09_run(cmp, g)
        expect_identical(out$components, cmp@components)
        expect_identical(length(out$edges), length(cmp@edges))
        for (side in 1:2) {
            source <- cellSummary(cmp, segmentationNames(cmp)[side])
            for (nm in colnames(source)) {
                expect_identical(out$cells[[side]][[nm]], source[[nm]])
            }
        }
        orig <- S4Vectors::mcols(cmp@edges)
        for (nm in colnames(orig)) {
            expect_identical(S4Vectors::mcols(out$edges)[[nm]], orig[[nm]])
        }
        # Identical polygons give every cross-side pair positive area, even ties.
        expect_identical(nrow(out$pairs), nrow(mat) * ncol(mat))
    }
})

test_that("relabeling preserves metrics and unique geometry matching", {
    skip_if_not_installed("sf")
    results <- lapply(list(list(c("x", "y"), c("a", "b")),
                         list(c("opaque-2", "opaque-1"), c("q", "z"))),
        function(ids) {
            cmp <- wp04_compare(diag(1L, 2), ids1 = ids[[1]], ids2 = ids[[2]])
            a <- wp08_box(0, 1, 0, 1)
            b <- wp08_box(2, 3, 0, 1)
            wp09_run(cmp, wp09_geometry(cmp, list(a, b), list(b, a)))
        })
    columns <- setdiff(colnames(results[[1]]$pairs), c("cell_id_1", "cell_id_2"))
    expect_identical(results[[1]]$pairs[, columns], results[[2]]$pairs[, columns])
    expect_identical(results[[1]]$summary, results[[2]]$summary)
})

test_that("no repair transformation or window operation is called", {
    skip_if_not_installed("sf")
    forbidden <- function(...) stop("Forbidden geometry mutation")
    testthat::local_mocked_bindings(st_make_valid = forbidden, st_buffer = forbidden,
        st_simplify = forbidden, st_transform = forbidden, st_crop = forbidden,
        st_set_crs = forbidden, st_union = forbidden, .package = "sf")
    cmp <- wp04_compare(matrix(1L, 1, 1))
    for (shape in list(wp08_box(), sf::st_polygon(),
                       sf::st_polygon(list(rbind(c(0, 0), c(1, 1), c(2, 2), c(0, 0)))))) {
        g <- wp09_geometry(cmp, list(shape), list(wp08_box()))
        before <- serialize(list(cmp, g), NULL)
        wp09_run(cmp, g)
        expect_identical(serialize(list(cmp, g), NULL), before)
    }
})

test_that("metric failures retain candidate rows and cannot certify matching", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(1L, 1, 1))
    g <- wp08_geometry(cmp)
    testthat::local_mocked_bindings(
        .geometryMeasurePair = function(...) stop("Backend metric failure"),
        .package = "SegStability")
    out <- wp09_run(cmp, g)
    expect_identical(out$pairs$geometry_status, "unchecked")
    expect_true(is.na(out$pairs$iou))
    expect_true(is.na(out$pairs$geometry_matched))
    expect_identical(out$summary$matching_status, "unchecked")
    expect_identical(out$validation$matching$status, "UNKNOWN")
    expect_true(is.na(out$summary$matching_objective))
})
