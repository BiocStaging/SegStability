# Assisted-by: OpenAI Codex
test_that("partition stores the frozen structure and exposes semantic accessors", {
    p <- wp02_partition()
    expect_s4_class(p, "SegmentationPartition")
    expect_true(methods::validObject(p))
    expect_identical(sort(methods::slotNames(p)), sort(c("name", "cells",
        "status", "cell_index", "ambiguity", "run_metadata", "provenance",
        "universe_fingerprint", "assignment_fingerprint", "fingerprint_version",
        "fingerprint_algorithm", "schema_version", "metadata")))
    expect_identical(cellCensus(p)$cell_id, c("c1", "c0", "unused"))
    expect_identical(assignmentStatus(p),
                     c("assigned", "unassigned", "ambiguous"))
    expect_identical(p@cell_index, c(1L, NA_integer_, NA_integer_))
    expect_identical(ambiguityData(p)$evidence$molecule_index, 3L)
    expect_identical(ambiguityData(p)$evidence$candidate_count, NA_integer_)
    expect_false(ambiguityData(p)$evidence$candidate_ids_known)
    expect_s4_class(runMetadata(p), "SimpleList")
    expect_s4_class(provenance(p), "SimpleList")
    expect_identical(p@schema_version, "1")
    expect_output(methods::show(p), "native")
    expect_false(any(c("molecule_id", "molecules", "geometry") %in%
                     methods::slotNames(p)))
    copied <- cellCensus(p)
    copied$cell_id[1L] <- "changed"
    expect_identical(cellCensus(p)$cell_id[1L], "c1")
})

test_that("partition constructor aligns external assignments exactly once", {
    p <- wp02_partition()
    q <- wp02_partition(molecule_id = c("m3", "m1", "m2"),
        status = c("ambiguous", "assigned", "unassigned"),
        cell_id = c(NA_character_, "c1", NA_character_))
    expect_identical(assignmentStatus(q), assignmentStatus(p))
    expect_identical(q@cell_index, p@cell_index)
    expect_identical(q@universe_fingerprint, p@universe_fingerprint)
    expect_identical(q@assignment_fingerprint, p@assignment_fingerprint)
    expect_error(wp02_partition(molecule_id = c("m1", "m1", "m3")),
                 "Duplicate")
    expect_error(wp02_partition(molecule_id = c("m1", "m2")), "missing")
    expect_error(wp02_partition(molecule_id = c("m1", "m2", "other")),
                 "missing.*extra")
    expect_error(wp02_partition(molecule_id = c("m1", "m2", "m3", "m4")),
                 "extra")
    expect_error(wp02_partition(molecule_id = c(1, 2, 3)), "unsupported")
    expect_error(wp02_partition(molecule_id = c("m1", NA_character_, "m3")),
                 "Missing")
    expect_error(wp02_partition(molecule_id = 1:3), "famil")
})

test_that("all certified molecule and cell families survive partition construction", {
    for (kind in c("integer", "character", "LLint", "integer64")) {
        if (kind == "integer64" && !requireNamespace("bit64", quietly = TRUE)) {
            next
        }
        ids <- id_lab_fixture(kind)$x
        n <- length(ids)
        p <- SegmentationPartition(
            molecules = S4Vectors::DataFrame(molecule_id = ids,
                feature_id = rep("g", n)), name = kind,
            cells = S4Vectors::DataFrame(cell_id = ids),
            molecule_id = ids[n:1], status = rep("assigned", n),
            cell_id = ids[n:1])
        expect_identical(p@cell_index, seq_len(n))
        expect_id_lab_exact(cellCensus(p)$cell_id, ids)
        path <- tempfile()
        saveRDS(p, path)
        q <- readRDS(path)
        unlink(path)
        expect_id_lab_exact(cellCensus(q)$cell_id, ids)
        expect_identical(q@assignment_fingerprint, p@assignment_fingerprint)
        expect_silent(.checkPartitionIntegrity(q))
    }
})

test_that("constructor rejects census and cross-field assignment failures", {
    expect_error(wp02_partition(cells = S4Vectors::DataFrame(
        cell_id = c("c1", "c1"))), "Duplicate")
    expect_error(wp02_partition(cells = S4Vectors::DataFrame(
        cell_id = c("c1", NA_character_))), "Missing")
    expect_error(wp02_partition(cells = S4Vectors::DataFrame(cell_id = 1.0)),
                 "unsupported")
    expect_error(wp02_partition(cell_id = c("other", NA, NA)), "Unknown")
    expect_error(wp02_partition(cell_id = c(NA, NA, NA)), "assigned")
    expect_error(wp02_partition(cell_id = c("c1", "c0", NA)), "forced")
    expect_error(wp02_partition(cell_id = c("c1", NA, "c0")), "forced")
    expect_error(wp02_partition(status = c("assigned", "unknown", "ambiguous")),
                 "status")
    expect_error(wp02_partition(status = "assigned"), "length")
    expect_error(wp02_partition(name = ""), "name")
    expect_error(wp02_partition(name = NA_character_), "name")
})

test_that("known candidates, counts and unknown ambiguity remain distinct", {
    p <- wp02_partition(ambiguity = wp02_ambiguity())
    a <- ambiguityData(p)
    expect_identical(a$evidence$candidate_count, 2L)
    expect_true(a$evidence$candidate_ids_known)
    expect_identical(S4Vectors::from(a$candidates), c(3L, 3L))
    expect_identical(S4Vectors::to(a$candidates), 1:2)
    b <- wp02_ambiguity()
    b$candidates <- b$candidates[2:1, ]
    expect_identical(wp02_partition(ambiguity = b)@assignment_fingerprint,
                     p@assignment_fingerprint)
    b <- wp02_ambiguity()
    b$candidates <- b$candidates[integer(), ]
    b$evidence$candidate_ids_known <- FALSE
    count_only <- wp02_partition(ambiguity = b)
    expect_identical(ambiguityData(count_only)$evidence$candidate_count, 2L)
    expect_false(ambiguityData(count_only)$evidence$candidate_ids_known)
    expect_identical(length(ambiguityData(count_only)$candidates), 0L)
    expect_false(identical(count_only@assignment_fingerprint,
                           wp02_partition()@assignment_fingerprint))
})

test_that("historical empty census and unresolved population cases are retained", {
    # Identity/status subset of oracle 144, not its graph calculations.
    m <- S4Vectors::DataFrame(molecule_id = integer(), feature_id = character())
    empty <- SegmentationPartition(m, "empty",
        S4Vectors::DataFrame(cell_id = integer()), integer(), character())
    expect_true(methods::validObject(empty))
    expect_identical(assignmentStatus(empty), character())
    m <- S4Vectors::DataFrame(molecule_id = 1:2, feature_id = c("g", "g"))
    zero <- SegmentationPartition(m, "zero",
        S4Vectors::DataFrame(cell_id = integer()), 1:2,
        c("unassigned", "ambiguous"))
    expect_identical(assignmentStatus(zero), c("unassigned", "ambiguous"))
    expect_identical(nrow(cellCensus(zero)), 0L)
    expect_identical(S4Vectors::nRnode(ambiguityData(zero)$candidates), 0L)
    # Final pre-API audit: cell 0 is a valid local ID, not a global sentinel.
    p <- SegmentationPartition(m, "zero-id",
        S4Vectors::DataFrame(cell_id = c(0L, -1L, 5L)), 1:2,
        c("assigned", "assigned"), c(0L, -1L))
    expect_identical(p@cell_index, 1:2)
    expect_identical(cellCensus(p)$cell_id, c(0L, -1L, 5L))
})

test_that("structural validity and explicit internal integrity checks are separate", {
    p <- wp02_partition(ambiguity = wp02_ambiguity())
    for (slot in c("name", "universe_fingerprint", "fingerprint_version",
                   "fingerprint_algorithm", "schema_version")) {
        bad <- p
        methods::slot(bad, slot) <- ""
        expect_error(methods::validObject(bad))
    }
    for (s in list(c(1L, 2L), c(1L, 4L, 3L), c(1L, NA_integer_, 3L))) {
        bad <- p
        bad@status <- s
        expect_error(methods::validObject(bad))
    }
    for (i in list(c(99L, NA_integer_, NA_integer_), c(1L, 2L, NA_integer_),
                   c(1L, NA_integer_, 2L))) {
        bad <- p
        bad@cell_index <- i
        expect_error(methods::validObject(bad))
    }
    bad <- p
    bad@ambiguity$candidates <- S4Vectors::Hits(1L, 1L, 3L, 3L)
    expect_error(methods::validObject(bad))
    bad <- p
    bad@ambiguity$evidence$candidate_count <- 1L
    expect_error(methods::validObject(bad))
    bad <- p
    bad@ambiguity$evidence$molecule_index <- 4L
    expect_error(methods::validObject(bad))
    bad <- p
    bad@ambiguity$evidence$candidate_ids_known <- NA
    expect_error(methods::validObject(bad))
    bad <- p
    bad@ambiguity$candidates <- S4Vectors::Hits(3L, 4L, 3L, 4L)
    expect_error(methods::validObject(bad))
    bad <- p
    bad@ambiguity$candidates <- S4Vectors::Hits(c(3L, 3L), c(1L, 1L), 3L, 3L)
    expect_error(methods::validObject(bad))
    bad <- p
    bad@cell_index[1L] <- 2L
    expect_true(methods::validObject(bad))
    expect_error(.checkPartitionIntegrity(bad), "assignment fingerprint")
    u <- .prepareUniverse(wp02_molecules())
    expect_silent(.checkPartitionIntegrity(p, u))
    u$molecules$molecule_id[1L] <- "other"
    u$fingerprint_value <- .universeFingerprint(u$molecules,
                                               u$features)$fingerprint_value
    expect_error(.checkPartitionIntegrity(p, u), "universe")
})

test_that("provenance unknowns and supplied metadata stay explicit", {
    p <- wp02_partition(provenance = S4Vectors::SimpleList(source = "fixture"),
        run_metadata = S4Vectors::SimpleList(seed = 7L),
        metadata = S4Vectors::SimpleList(note = "retained"))
    expect_identical(provenance(p)$source, "fixture")
    expect_identical(provenance(p)$software_version, NA_character_)
    expect_identical(runMetadata(p)$seed, 7L)
    expect_identical(p@metadata$note, "retained")
    expect_error(wp02_partition(provenance = list()), "SimpleList")
    expect_error(wp02_partition(run_metadata = list()), "SimpleList")
})
