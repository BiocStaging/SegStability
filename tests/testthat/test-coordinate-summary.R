# Assisted-by: OpenAI Codex
test_that("summaries expose exact stored coordinate declarations without gates", {
    records <- list(S4Vectors::SimpleList(),
        S4Vectors::SimpleList(units = "micrometres", origin = NA_character_,
            registration = "UNKNOWN", unknown_integer = NA_integer_,
            unknown_real = NA_real_, unknown_logical = NA),
        S4Vectors::SimpleList(units = "micrometres", origin = c(0, 0),
            axis_directions = c("right", "up"),
            transforms = S4Vectors::SimpleList(matrix = diag(3),
                history = list(list(name = "declared", parameters = c(1L, 2L)))),
            registration = "declared identity", specimen_id = "specimen"))
    for (info in records) {
        x <- wp03_set(coordinate_info = info)
        cmp <- compareSegmentations(x, "native", "other")
        before <- serialize(list(x, cmp), NULL)
        expect_identical(summary(x)$coordinate_info, info)
        expect_identical(summary(cmp)$coordinate_info, info)
        f <- tempfile()
        saveRDS(list(x, cmp), f)
        restored <- readRDS(f)
        unlink(f)
        expect_identical(summary(restored[[1]])$coordinate_info, info)
        expect_identical(summary(restored[[2]])$coordinate_info, info)
        returned <- summary(x)$coordinate_info
        returned$units <- "changed"
        expect_identical(serialize(list(x, cmp), NULL), before)
    }
    testthat::local_mocked_bindings(
        validateCoordinates = function(...) stop("coordinate gate invoked"),
        validateGeometry = function(...) stop("geometry gate invoked"),
        .checkSetIntegrity = function(...) stop("integrity gate invoked"))
    expect_identical(summary(x)$coordinate_info, info)
    expect_identical(summary(cmp)$coordinate_info, info)
})
