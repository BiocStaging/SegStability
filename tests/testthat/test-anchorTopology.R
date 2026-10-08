test_that("anchorTopology has exactly the frozen public signature", {
    expect_identical(names(formals(anchorTopology)),
                     c("comparison", "geometry", "window", "coordinate_info"))
    expect_null(formals(anchorTopology)$coordinate_info)
    expect_false("asSpatialExperiment" %in% getNamespaceExports("SegStability"))
})

test_that("topology validates the exact rectangle without bounds repair", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare()
    g <- wp08_geometry(cmp)
    bad <- list(NULL, c(0, 10, 0, 10), c(xmin = 0, xmax = 10, ymin = 0),
        c(xmin = 0, xmax = 10, ymin = 0, ymin = 10),
        c(xmin = 10, xmax = 0, ymin = 0, ymax = 10),
        c(xmin = 0, xmax = 10, ymin = 1, ymax = 1),
        c(xmin = 0, xmax = Inf, ymin = 0, ymax = 10),
        c(xmin = 0, xmax = NA, ymin = 0, ymax = 10),
        structure(wp08_window(), class = "bbox"),
        structure(wp08_window(), dim = c(2L, 2L)), as.list(wp08_window()))
    for (w in bad) expect_error(anchorTopology(cmp, g, w, wp08_info()), "window")
    expect_error(anchorTopology(cmp, g, coordinate_info = wp08_info()), "window")
    w <- wp08_window()[c(4, 1, 3, 2)]
    expect_identical(anchorTopology(cmp, g, w, wp08_info())$window, w)
    expect_identical(anchorTopology(cmp, g, as.integer(wp08_window()) |>
        setNames(names(wp08_window())), wp08_info())$anchors$anchor_topology,
        c("one_to_one_like", "one_to_one_like"))
})

test_that("output preserves complete censuses graphs ordering and serialization", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(c(3, 0, 0, 0), 2))
    g <- wp08_geometry(cmp)
    r <- anchorTopology(cmp, g, wp08_window(), wp08_info())
    expect_s4_class(r, "SimpleList")
    expect_identical(names(r), c("anchors", "cells", "edges", "validation", "window"))
    expect_identical(names(r$cells), segmentationNames(cmp))
    expect_s4_class(r$anchors, "DataFrame")
    expect_s4_class(r$edges, "Hits")
    for (i in 1:2) {
        original <- cellSummary(cmp, segmentationNames(cmp)[i])
        expect_identical(r$cells[[i]][, colnames(original)], original)
        expect_identical(r$cells[[i]]$selection_stratum,
                         c("interior_anchor", "zero_correspondence"))
    }
    expect_identical(S4Vectors::from(r$edges), S4Vectors::from(correspondenceEdges(cmp)))
    expect_identical(S4Vectors::nLnode(r$edges), 2L)
    shuffled <- lapply(g, function(x) x[2:1, ])
    expect_identical(anchorTopology(cmp, shuffled, wp08_window(), wp08_info())$anchors,
                     r$anchors)
    expect_identical(anchorTopology(cmp, shuffled, wp08_window(), wp08_info())$cells, r$cells)
    file <- tempfile(fileext = ".rds")
    on.exit(unlink(file))
    saveRDS(r, file)
    expect_identical(readRDS(file), r)
    wp08_unchanged(cmp, g, anchorTopology(cmp, g, wp08_window(), wp08_info()))
    for (m in list(matrix(0, 0, 0), matrix(0, 2, 0), matrix(0, 2, 2))) {
        empty <- wp08_run(m)
        expect_identical(nrow(empty$anchors), 0L)
        expect_identical(length(empty$edges), 0L)
        expect_identical(nrow(empty$cells[[1]]), nrow(m))
        expect_identical(nrow(empty$cells[[2]]), ncol(m))
        expect_type(empty$anchors$anchor_topology, "character")
    }
})

test_that("spatial analysis requires declarations and optional sf explicitly", {
    cmp <- wp04_compare()
    local_mocked_bindings(.geometryAvailable = function() FALSE)
    expect_error(anchorTopology(cmp, list(left = NULL, right = NULL),
                               wp08_window(), wp08_info()), "sf")
})

test_that("known common declarations are required without source authenticity claims", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare()
    g <- wp08_geometry(cmp)
    expect_error(anchorTopology(cmp, g, wp08_window()), "prerequisite")
    for (field in names(wp08_info())) {
        info <- wp08_info()
        info[[field]] <- NULL
        expect_error(anchorTopology(cmp, g, wp08_window(), info), "prerequisite")
    }
    info <- wp08_info()
    info$specimen_id <- "other"
    expect_error(anchorTopology(cmp, g, wp08_window(), info), "prerequisite")
    cmp@coordinate_info <- wp08_info()
    expect_identical(anchorTopology(cmp, g, wp08_window())$anchors$anchor_topology,
                     c("one_to_one_like", "one_to_one_like"))
    info <- wp08_info()
    info$units <- "pixels"
    expect_error(anchorTopology(cmp, g, wp08_window(), info), "prerequisite")
    r <- anchorTopology(cmp, g, wp08_window(), wp08_info())
    expect_identical(validation_row(r$validation$coordinates,
        "topology.registration_evidence")$status, "UNKNOWN")
    bad <- cmp
    S4Vectors::mcols(bad@edges)$unique_dominant_1 <- rep(FALSE, length(bad@edges))
    expect_error(anchorTopology(bad, g, wp08_window(), wp08_info()), "prerequisite")
})
