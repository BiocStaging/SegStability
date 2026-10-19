# Assisted-by: OpenAI Codex
test_that("synthetic example preserves its documented census and masses", {
    data("segmentationExample", package = "SegStability")
    d <- segmentationExample
    expect_identical(nrow(d$molecules), 24L)
    expect_identical(length(unique(d$molecules$molecule_id)), 24L)
    expect_identical(vapply(d$cells, nrow, integer(1)),
                     c(first = 6L, second = 7L))
    expect_identical(d$ambiguity$first$evidence$candidate_count, 2L)

    make_partition <- function(side) {
        a <- d$assignments[[side]]
        ambiguity <- d$ambiguity[[side]]
        if (!is.null(ambiguity)) {
            ambiguity <- S4Vectors::SimpleList(
                evidence = S4Vectors::DataFrame(ambiguity$evidence),
                candidates = S4Vectors::DataFrame(ambiguity$candidates))
        }
        SegmentationPartition(d$molecules, side, d$cells[[side]],
                              a$molecule_id, a$status, a$cell_id,
                              ambiguity = ambiguity)
    }
    ss <- SegmentationSet(d$molecules,
                          list(first = make_partition("first"),
                               second = make_partition("second")),
                          specimen_id = "synthetic-example")
    cmp <- compareSegmentations(ss, "first", "second")
    expect_identical(sum(S4Vectors::mcols(correspondenceEdges(cmp))$weight),
                     20)
    expect_identical(sum(cellSummary(cmp, "first")$assigned_mass), 21)
    expect_identical(sum(cellSummary(cmp, "second")$assigned_mass), 22)
    expect_setequal(as.character(componentSummary(cmp)$component_class),
                    c("one_to_one", "split", "merge", "unmatched"))
    states <- stateStability(cmp, d$states)
    expect_identical(states$global$resolved_mass, 18)
    expect_identical(states$global$same_state_mass, 15)
})
