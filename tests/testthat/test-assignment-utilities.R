test_that("assignment utilities retain exact status and cell semantics", {
    s <- .canonicalStatus(c("assigned", "unassigned", "ambiguous"))
    expect_identical(s, 1:3)
    expect_error(.canonicalStatus(c("assigned", NA_character_)), "status")
    expect_error(.canonicalStatus("unknown"), "status")
    expect_error(.canonicalStatus(1:3), "status")
    expect_identical(.mapCells(c("zero", NA, NA), c("zero", "unused"), s),
                     c(1L, NA_integer_, NA_integer_))
    expect_error(.mapCells(c(NA, NA, NA), "zero", s), "assigned")
    expect_error(.mapCells(c("zero", "zero", NA), "zero", s), "forced")
    expect_error(.mapCells(c("other", NA, NA), "zero", s), "Unknown")
    expect_error(.mapCells(c(0L, NA_integer_, NA_integer_), "0", s), "famil")
})

test_that("ambiguity rejects malformed evidence without repairs", {
    make <- function() S4Vectors::SimpleList(
        evidence = S4Vectors::DataFrame(molecule_id = "m",
            candidate_count = 2L, candidate_ids_known = TRUE),
        candidates = S4Vectors::DataFrame(molecule_id = c("m", "m"),
            cell_id = c("c1", "c2")))
    check <- function(a) .prepareAmbiguity(a, "m", c("c1", "c2"), 3L)
    a <- make()
    a$candidates$cell_id[2L] <- "c1"
    expect_error(check(a), "Duplicate")
    a <- make()
    a$evidence$candidate_ids_known <- FALSE
    expect_error(check(a), "unknown")
    a <- make()
    a$evidence$molecule_id <- "other"
    expect_error(check(a), "ambiguous")
    for (count in list(-1L, 1.5, Inf, NaN, "2", factor("2"), TRUE)) {
        expect_error(.candidateCounts(count), "candidate_count")
        a <- make()
        a$evidence$candidate_count <- count
        expect_error(check(a), "candidate_count")
    }
})

test_that("ambiguity maps external IDs into full universe and census dimensions", {
    a <- S4Vectors::SimpleList(
        evidence = S4Vectors::DataFrame(molecule_id = c("m3", "m2"),
            candidate_count = c(3L, 2L), candidate_ids_known = c(FALSE, TRUE)),
        candidates = S4Vectors::DataFrame(molecule_id = c("m2", "m2"),
            cell_id = c("c3", "c1")))
    b <- .prepareAmbiguity(a, c("m1", "m2", "m3", "m4"),
                           c("c1", "c2", "c3"), c(1L, 3L, 3L, 3L))
    expect_identical(b$evidence$molecule_index, 2:4)
    expect_identical(b$evidence$candidate_count, c(2L, 3L, NA_integer_))
    expect_identical(b$evidence$candidate_ids_known, c(TRUE, FALSE, FALSE))
    expect_identical(S4Vectors::from(b$candidates), c(2L, 2L))
    expect_identical(S4Vectors::to(b$candidates), c(1L, 3L))
    expect_identical(S4Vectors::nLnode(b$candidates), 4L)
    expect_identical(S4Vectors::nRnode(b$candidates), 3L)
    bad <- a
    bad$evidence$candidate_count[2L] <- 1L
    expect_error(.prepareAmbiguity(bad, c("m1", "m2", "m3", "m4"),
        c("c1", "c2", "c3"), c(1L, 3L, 3L, 3L)), "count")
    bad <- a
    bad$candidates$cell_id[1L] <- "unknown"
    expect_error(.prepareAmbiguity(bad, c("m1", "m2", "m3", "m4"),
        c("c1", "c2", "c3"), c(1L, 3L, 3L, 3L)), "Unknown")
    bad <- a
    bad$candidates$molecule_id[1L] <- "m1"
    expect_error(.prepareAmbiguity(bad, c("m1", "m2", "m3", "m4"),
        c("c1", "c2", "c3"), c(1L, 3L, 3L, 3L)), "ambiguous")
})
