# Assisted-by: OpenAI Codex
test_that("universe fingerprints are deterministic, typed and order-sensitive", {
    make <- function(ids, features = rep("g", length(ids))) {
        .prepareUniverse(S4Vectors::DataFrame(
            molecule_id = ids, feature_id = features))
    }
    a <- make(c("a", "b"))
    expect_identical(a$fingerprint_version, "1")
    expect_identical(a$fingerprint_algorithm, "sha256")
    expect_match(a$fingerprint_value, "^[0-9a-f]{64}$")
    expect_identical(a$fingerprint_value, make(c("a", "b"))$fingerprint_value)
    for (b in list(make(c("b", "a")), make(c("a", "c")),
                   make(c("a", "b"), c("g", "h")),
                   make(c("a", "b"), c("h", "g")))) {
        expect_false(identical(a$fingerprint_value, b$fingerprint_value))
    }
    expect_false(identical(make(1:2)$fingerprint_value,
                           make(c("1", "2"))$fingerprint_value))
    expect_false(identical(make(c("a|b", "c"))$fingerprint_value,
                           make(c("a", "b|c"))$fingerprint_value))
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(a, path)
    b <- readRDS(path)
    expect_identical(b$fingerprint_value, a$fingerprint_value)
    expect_silent(.validateUniverse(b))
    b$molecules$molecule_id[1L] <- "corrupted"
    expect_error(.validateUniverse(b), "fingerprint")
    b <- a
    b$fingerprint_version <- "unsupported"
    expect_error(.validateUniverse(b), "version")
})

test_that("canonical empty fields and missing values remain unambiguous", {
    expect_identical(.frameValues(character()), "")
    expect_false(identical(.frameValues(NA_character_), .frameValues("NA")))
    expect_false(identical(.frameValues(c("a", "b")), .frameValues("ab")))
})

test_that("64-bit universe fingerprints never use ordinary numeric conversion", {
    local_traces <- function() {
        trace(".from_LLint_to_numeric", where = asNamespace("S4Vectors"),
              tracer = quote(stop("unsafe double conversion")), print = FALSE)
        function() untrace(".from_LLint_to_numeric",
                           where = asNamespace("S4Vectors"))
    }
    cleanup <- local_traces()
    on.exit(cleanup())
    ids <- S4Vectors::as.LLint(id_lab_large_values())
    a <- .prepareUniverse(S4Vectors::DataFrame(
        molecule_id = ids, feature_id = rep("g", length(ids))))
    expect_silent(.validateUniverse(a))
    expect_identical(.alignMolecules(ids[9:1], ids), 9:1)
    expect_identical(.idOrder(ids), id_lab_fixture("LLint")$sorted)
    p <- SegmentationPartition(S4Vectors::DataFrame(
        molecule_id = ids, feature_id = rep("g", length(ids))), "LLint",
        S4Vectors::DataFrame(cell_id = ids), ids[9:1],
        rep("assigned", 9L), ids[9:1])
    expect_silent(.checkPartitionIntegrity(p, a))
    if (requireNamespace("bit64", quietly = TRUE)) {
        trace("as.double.integer64", where = asNamespace("bit64"),
              tracer = quote(stop("unsafe double conversion")), print = FALSE)
        on.exit(untrace("as.double.integer64", where = asNamespace("bit64")),
                add = TRUE)
        ids <- bit64::as.integer64(id_lab_large_values())
        b <- .prepareUniverse(S4Vectors::DataFrame(
            molecule_id = ids, feature_id = rep("g", length(ids))))
        expect_silent(.validateUniverse(b))
        expect_identical(.alignMolecules(ids[9:1], ids), 9:1)
        expect_identical(.idOrder(ids), id_lab_fixture("integer64")$sorted)
        expect_false(identical(a$fingerprint_value, b$fingerprint_value))
        p <- SegmentationPartition(S4Vectors::DataFrame(
            molecule_id = ids, feature_id = rep("g", length(ids))), "integer64",
            S4Vectors::DataFrame(cell_id = ids), ids[9:1],
            rep("assigned", 9L), ids[9:1])
        expect_silent(.checkPartitionIntegrity(p, b))
    }
})

test_that("assignment fingerprint binds universe, census, status and ambiguity", {
    cells <- S4Vectors::DataFrame(cell_id = c("c", "unused"))
    a <- .prepareAmbiguity(NULL, c("m1", "m2"), cells$cell_id, c(1L, 3L))
    fp <- function(u = "universe", c = cells, s = c(1L, 3L),
                   i = c(1L, NA_integer_), e = a) {
        .assignmentFingerprint(u, c, s, i, e)
    }
    original <- fp()
    expect_identical(original, fp())
    modified <- a
    modified$evidence$candidate_count <- 2L
    changed <- cells
    changed$cell_id[2L] <- "different-unused"
    for (b in list(fp(u = "other"), fp(c = changed), fp(e = modified),
                   fp(s = c(1L, 2L)), fp(i = c(2L, NA_integer_)))) {
        expect_false(identical(original$fingerprint_value, b$fingerprint_value))
    }
    expect_identical(original$fingerprint_version, "1")
    expect_identical(original$fingerprint_algorithm, "sha256")
})
