test_that("the geometry API and exact overlap quantities are explicit", {
    skip_if_not_installed("sf")
    expect_identical(names(formals(compareGeometry)),
                     c("comparison", "geometry", "coordinate_info"))
    cmp <- wp04_compare(matrix(3L, 1, 1))
    box <- wp08_box(0, 2, 0, 2)
    others <- list(box, wp08_box(3, 5, 0, 2), wp08_box(1, 3, 0, 2),
                   wp08_box(0, 1, 0, 1))
    expected <- list(c(4, 4, 4, 4, 1, 1, 1, 1),
        c(4, 4, 0, 8, 0, 0, 0, 0), c(4, 4, 2, 6, 1/3, .5, .5, .5),
        c(4, 1, 1, 4, .25, .4, .25, 1))
    metrics <- c("area_1", "area_2", "intersection_area", "union_area",
                 "iou", "dice", "overlap_fraction_1", "overlap_fraction_2")
    for (i in seq_along(others)) {
        g <- wp09_geometry(cmp, list(box), list(others[[i]]))
        before <- serialize(list(cmp, g), NULL)
        result <- wp09_run(cmp, g)
        expect_identical(names(result), c("pairs", "matching", "cells", "edges",
            "components", "summary", "validation", "coordinate_info"))
        expect_equal(as.numeric(unlist(as.list(result$pairs[1, metrics]))), expected[[i]])
        expect_identical(serialize(list(cmp, g), NULL), before)
        expect_identical(result$components, componentSummary(cmp))
    }
})

test_that("geometry-only matches cross original components without changing them", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(diag(2L, 2))
    a <- wp08_box(0, 1, 0, 1)
    b <- wp08_box(2, 3, 0, 1)
    g <- wp09_geometry(cmp, list(a, b), list(b, a))
    before <- serialize(cmp, NULL)
    out <- wp09_run(cmp, g)
    expect_identical(nrow(out$pairs), 4L)
    expect_identical(out$matching$edge_index, c(NA_integer_, NA_integer_))
    expect_true(all(out$matching$component_id_1 != out$matching$component_id_2))
    expect_identical(out$matching$cell_index_1, 1:2)
    expect_identical(out$matching$cell_index_2, 2:1)
    expect_identical(out$components, cmp@components)
    expect_identical(serialize(cmp, NULL), before)
    original <- S4Vectors::mcols(cmp@edges)
    annotated <- S4Vectors::mcols(out$edges)
    for (nm in colnames(original)) expect_identical(annotated[[nm]], original[[nm]])
    expect_identical(S4Vectors::from(out$edges), S4Vectors::from(cmp@edges))
    expect_identical(S4Vectors::to(out$edges), S4Vectors::to(cmp@edges))
    expect_identical(S4Vectors::nLnode(out$edges), S4Vectors::nLnode(cmp@edges))
    expect_identical(S4Vectors::nRnode(out$edges), S4Vectors::nRnode(cmp@edges))
})
