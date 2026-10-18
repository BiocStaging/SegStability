test_that("the frozen class and internal constructor represent census isolates", {
    args <- wp04_class_input()
    x <- do.call(SegStability:::.newComparison, args)
    expect_s4_class(x, "SegmentationComparison")
    expect_true(validObject(x))
    expect_identical(sort(slotNames(x)), sort(names(args)))
    expect_identical(segmentationNames(x), c("left", "right"))
    expect_identical(correspondenceEdges(x), args$edges)
    expect_identical(cellSummary(x, "left"), args$cell_summary_1)
    expect_identical(componentSummary(x), args$components)
    expect_identical(assignmentSummary(x), args$assignment_summary)
    expect_error(cellSummary(x), "segmentation")
    expect_error(cellSummary(x, "missing"), "segmentation")
    expect_error(cellSummary(x, c("left", "right")), "segmentation")
    expect_identical(names(summary(x)), c("segmentation_names", "universe_size",
        "cell_summary_1", "cell_summary_2", "components", "assignment_summary",
        "coordinate_info"))
    expect_output(show(x), "SegmentationComparison")
    file <- tempfile(fileext = ".rds")
    on.exit(unlink(file))
    saveRDS(x, file)
    expect_identical(readRDS(file), x)
    y <- cellSummary(x, "left")
    y$assigned_mass[1L] <- 99
    expect_identical(cellSummary(x, "left")$assigned_mass, c(2, 0))
})

test_that("malformed direct comparison objects fail structural validity", {
    x <- do.call(SegStability:::.newComparison, wp04_class_input())
    bad <- function(slot, value) {
        y <- x
        slot(y, slot) <- value
        expect_false(isTRUE(validObject(y, test = TRUE)))
    }
    bad("segmentation_names", c("left", "left"))
    bad("schema_version", "2")
    bad("assignment_fingerprints", "bad")
    bad("universe_size", 2^53 + 2)
    bad("universe_size", -1)
    h <- x@edges
    for (w in c(0, -1, 0.5, Inf, NA_real_, 2^53 + 2)) {
        S4Vectors::mcols(h)$weight <- w
        bad("edges", h)
    }
    h <- x@edges[c(1L, 1L)]
    bad("edges", h)
    h <- S4Vectors::Hits(1L, 1L, 3L, 2L)
    S4Vectors::mcols(h) <- S4Vectors::mcols(x@edges)
    bad("edges", h)
    c <- x@cell_summary_1
    c$matched_mass[1L] <- 3
    bad("cell_summary_1", c)
    c <- x@cell_summary_1
    c$component_id[1L] <- 9L
    bad("cell_summary_1", c)
    c <- x@components
    c$component_class[1L] <- "one_to_one_like"
    bad("components", c)
    a <- x@assignment_summary
    a$n_molecules[1L] <- 3
    bad("assignment_summary", a)
    for (column in c("fraction_1", "dominant_1", "reciprocal_dominant")) {
        h <- x@edges
        S4Vectors::mcols(h)[[column]] <- NA
        bad("edges", h)
    }
    h <- x@edges
    S4Vectors::mcols(h)$fraction_1 <- 0.5
    bad("edges", h)
    y <- wp04_compare()
    y@edges <- y@edges[2:1]
    expect_false(isTRUE(validObject(y, test = TRUE)))
})

test_that("exact counts support the frozen upper boundary without ID coercion", {
    args <- wp04_class_input()
    args$universe_size <- 2^53
    args$assignment_summary$n_molecules[1L] <- 2^53
    for (nm in c("cell_summary_1", "cell_summary_2")) {
        args[[nm]]$assigned_mass[1L] <- 2^53
        args[[nm]]$matched_mass[1L] <- 2^53
    }
    S4Vectors::mcols(args$edges)$weight <- 2^53
    expect_true(validObject(do.call(SegStability:::.newComparison, args)))
    args$cell_summary_1$unmatched_mass[1L] <- 1
    expect_error(do.call(SegStability:::.newComparison, args), "masses")
    args$cell_summary_1$unmatched_mass[1L] <- 0
    # Rounded sum(2^53, 1) must not hide a universe conservation violation.
    args$assignment_summary$n_molecules[9L] <- 1
    expect_error(do.call(SegStability:::.newComparison, args), "conservation")
    expect_error(SegStability:::.graphComponents(integer(), integer(),
        .Machine$integer.max, 1L), "positional index range")
})
