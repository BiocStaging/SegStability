test_that("production identity paths preserve all Gate D families exactly", {
    kinds <- c("integer", "character", "LLint", "integer64")
    for (kind in kinds) {
        if (kind == "integer64" && !requireNamespace("bit64", quietly = TRUE)) {
            next
        }
        f <- id_lab_fixture(kind)
        expect_id_lab_exact(.validateIds(f$x, "molecule_id"), f$x)
        expect_identical(.idOrder(f$x), f$sorted)
        expect_identical(.matchIds(f$x[rev(seq_along(f$x))], f$x),
                         rev(seq_along(f$x)))
        expect_error(.validateIds(f$x[c(1L, 1L)], "cell_id"), "Duplicate")
        expect_id_lab_exact(S4Vectors::DataFrame(id = f$x)$id, f$x)
    }
})

test_that("production ID validation rejects unsupported paths", {
    invalid <- list(NULL, 1, TRUE, factor("a"), 1 + 1i, as.raw(1L),
                    NA_integer_, NA_character_, NaN, Inf, "",
                    structure(1L, class = "custom_id"), list(1L, "1"))
    for (x in invalid) {
        expect_error(.validateIds(x, "molecule_id"))
    }
    expect_error(.matchIds(1L, "1"), "famil")
    expect_error(.matchIds(S4Vectors::as.LLint("1"), 1L), "famil")
    expect_identical(.matchIds(c("001", " 1 "), c("1", "001", " 1 ")),
                     c(2L, 3L))
})

test_that("external molecule alignment rejects mismatch and returns exact order", {
    u <- c("m1", "m2", "m3")
    expect_identical(.alignMolecules(u[c(3L, 1L, 2L)], u), c(2L, 3L, 1L))
    expect_error(.alignMolecules(c("m1", "m1", "m3"), u), "Duplicate")
    expect_error(.alignMolecules(u[-1L], u), "missing")
    expect_error(.alignMolecules(c(u, "m4"), u), "extra")
    expect_error(.alignMolecules(c("m1", "m2", "other"), u), "missing.*extra")
    expect_error(.alignMolecules(1:3, u), "famil")
})

test_that("feature IDs have their own character-only opaque contract", {
    ids <- c("g", " g", "g", "G", "g ", "\u00e9", "e\u0301")
    u <- .prepareUniverse(S4Vectors::DataFrame(
        molecule_id = seq_along(ids), feature_id = ids))
    expect_identical(u$features$feature_id, ids[c(1L, 2L, 4L, 5L, 6L, 7L)])
    expect_identical(u$molecules$feature_index, c(1L, 2L, 1L, 3L, 4L, 5L, 6L))
    expect_identical(u$molecules$molecule_id, seq_along(ids))
    expect_false("feature_id" %in% colnames(u$molecules))
    for (bad in list(1L, 1, TRUE, factor("g"), NA_character_, "",
                    S4Vectors::as.LLint("1"))) {
        expect_error(.prepareUniverse(S4Vectors::DataFrame(
            molecule_id = 1L, feature_id = bad)), "feature_id")
    }
    bad <- u
    bad$features$feature_id[2L] <- bad$features$feature_id[1L]
    expect_error(.validateUniverse(bad), "Duplicate")
    bad <- u
    bad$molecules$feature_index[1L] <- 99L
    expect_error(.validateUniverse(bad), "feature_index")
})

test_that("universe validation retains optional columns without ID copies", {
    m <- S4Vectors::DataFrame(molecule_id = c("m2", "m1"),
        feature_id = c("g", "g"), x = c(2.5, 1.5), kept = c(TRUE, FALSE))
    u <- .prepareUniverse(m)
    expect_identical(u$molecules$x, m$x)
    expect_identical(u$molecules$kept, m$kept)
    expect_identical(u$molecules$molecule_index, 1:2)
    expect_identical(u$molecules$molecule_id, m$molecule_id)
    expect_error(.prepareUniverse(S4Vectors::DataFrame(
        molecule_id = c("m", "m"), feature_id = c("g", "g"))), "Duplicate")
    expect_error(.prepareUniverse(S4Vectors::DataFrame(
        molecule_id = c("m", NA_character_), feature_id = c("g", "g"))),
        "Missing")
    expect_error(.prepareUniverse(S4Vectors::DataFrame(molecule_id = "m")),
                 "feature_id")
})
