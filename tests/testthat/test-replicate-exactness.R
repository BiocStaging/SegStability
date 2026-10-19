# Assisted-by: OpenAI Codex
test_that("run-local labels census order and molecule input order do not change results", {
    reference <- replicateAgreement(wp10_set(), c("run_a", "run_b"))
    for (kind in c("integer", "character", "LLint", "integer64")) {
        if (kind == "integer64") skip_if_not_installed("bit64")
        ids <- id_lab_fixture(kind)$x[1:3]
        x <- wp10_set(cell_ids = list(ids, rev(ids)),
            census_orders = list(c(3L, 1L, 2L), c(2L, 3L, 1L)),
            molecule_orders = list(c(4L, 2L, 1L, 3L), c(2L, 4L, 3L, 1L)))
        expect_identical(replicateAgreement(x, c("run_a", "run_b")), reference)
        expect_id_lab_exact(cellCensus(x@segmentations$run_a)$cell_id, ids[c(3, 1, 2)])
        path <- tempfile()
        saveRDS(x, path)
        restored <- readRDS(path)
        unlink(path)
        expect_identical(replicateAgreement(restored, c("run_a", "run_b")), reference)
    }
    relabeled <- wp10_set(cell_ids = list(c("same", "other", "zero"),
                                         c("other", "same", "empty")))
    expect_identical(replicateAgreement(relabeled, c("run_a", "run_b")), reference)
})

test_that("molecule identity stays exact across all audited ID families", {
    for (kind in c("integer", "character", "LLint", "integer64")) {
        if (kind == "integer64") skip_if_not_installed("bit64")
        ids <- id_lab_fixture(kind)$x
        n <- length(ids)
        a <- rep(c(1L, 2L), length.out = n)
        x <- wp10_set(list(a = a, b = a), molecule_ids = ids,
                     molecule_orders = list(seq_len(n), rev(seq_len(n))))
        expect_id_lab_exact(moleculeData(x)$molecule_id, ids)
        expect_identical(replicateAgreement(x, c("a", "b"))$pairwise_summary$
            ari_jointly_assigned, 1)
    }
})

test_that("choose2 accumulation stays exact beyond double integer precision", {
    expect_identical(as.character(.replicateChoose2(.Machine$integer.max)),
                     "2305843005992468481")
    expect_identical(as.character(.replicateChoose2(c(1073741824L, 1073741823L))),
                     "1152921502459363329")
    expect_identical(as.character(.replicateChoose2(integer())), "0")
    expect_identical(as.character(.replicateChoose2(c(0L, 1L, 2L, 3L))), "4")
})

test_that("ARI zero denominator decisions precede floating conversion", {
    total <- wp10_ll("2305843005992468481")
    adjacent <- total - 1L
    expect_identical(suppressWarnings(as.numeric(total)), suppressWarnings(as.numeric(adjacent)))
    expect_identical(.replicateARIFromCounts(total, total, total, total), NA_real_)
    # Arithmetic-only adjacent terms force a nonzero denominator hidden by doubles.
    expect_identical(.replicateARIFromCounts(adjacent, adjacent, adjacent, total), 1)
    margin <- wp10_ll("1152921502459363329")
    expect_identical(.replicateARIFromCounts(margin, margin, margin, total), 1)
    expect_identical(.replicateARIFromCounts(margin, total, margin, total), 0)
    expect_identical(.replicateARIFromCounts(wp10_ll("0"), wp10_ll("0"),
                                             wp10_ll("0"), total), NA_real_)
})

test_that("large sparse contingencies produce deterministic numeric ARI", {
    h <- S4Vectors::Hits(c(1L, 2L), c(1L, 2L), 3L, 3L)
    S4Vectors::mcols(h) <- S4Vectors::DataFrame(weight = c(1073741824, 1073741823))
    expect_identical(.replicateARI(h), 1)
    # A one-molecule split of a near-limit one-cell partition has ARI exactly zero.
    h <- S4Vectors::Hits(c(1L, 1L), c(1L, 2L), 1L, 2L)
    S4Vectors::mcols(h) <- S4Vectors::DataFrame(weight = c(2147483646, 1))
    expect_identical(.replicateARI(h), 0)
    h <- S4Vectors::Hits(c(1L, 1L, 2L, 2L), c(1L, 2L, 1L, 2L), 2L, 2L)
    S4Vectors::mcols(h) <- S4Vectors::DataFrame(
        weight = c(600000000, 300000000, 200000000, 900000000))
    # Independent Python integer/Fraction calculation, with all pair sums >2^53:
    # exact ARI = 3327999992 / 13327999987.
    expected <- 3327999992 / 13327999987
    expect_equal(.replicateARI(h), expected, tolerance = 1e-15)
    expect_identical(.replicateARI(h), .replicateARI(h))
})
