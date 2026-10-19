test_that("saved sets expose exact named bound partitions", {
    x <- wp03_set()
    before <- serialize(x, NULL)
    p <- x[["native"]]
    expect_identical(p, x@segmentations$native)
    expect_identical(p@name, "native")
    expect_identical(cellCensus(p)$cell_id, c("c1", "c0", "unused"))
    expect_identical(assignmentStatus(p), c("assigned", "unassigned", "ambiguous"))
    expect_identical(provenance(p), provenance(x@segmentations$native))
    expect_identical(runMetadata(p), runMetadata(x@segmentations$native))
    expect_identical(ambiguityData(p), ambiguityData(x@segmentations$native))
    expect_identical(p@universe_fingerprint, x@universe_fingerprint)
    file <- tempfile()
    on.exit(unlink(file))
    saveRDS(x, file)
    expect_identical(readRDS(file)[["native"]], p)
    p@cells$cell_id[1L] <- "changed"
    expect_identical(serialize(x, NULL), before)
    expect_error(x[["native"]] <- p)
    expect_error(x["native"])
    expect_error(as.list(x))
})

test_that("partition retrieval rejects selectors outside exact names", {
    x <- wp03_set()
    for (selector in list("nat", "missing", "", character(), NA_character_,
                         c("native", "other"), 1L, 1, NULL, TRUE,
                         matrix("native"))) {
        expect_error(x[[selector]], "segmentation")
    }
    expect_error(x[[]], "segmentation")
    expect_error(x[["native", 1L]], "argument")
    expect_error(x[["native", exact = FALSE]], "exact")
    expect_error(x[["native", exact = NA]], "exact")
    expect_error(x[["native", extra = TRUE]], "argument")
    expect_identical(x[["native", exact = TRUE]], x[["native"]])
})

test_that("retrieval preserves all native census ID families and integrity", {
    for (kind in c("integer", "character", "LLint", "integer64")) {
        if (kind == "integer64" && !requireNamespace("bit64", quietly = TRUE)) next
        ids <- id_lab_fixture(kind)$x
        p <- wp02_partition(cells = S4Vectors::DataFrame(cell_id = ids),
            cell_id = ids[c(1L, 1L, 1L)], status = rep("assigned", 3L))
        x <- SegmentationSet(wp02_molecules(), list(native = p), "specimen")
        expect_id_lab_exact(cellCensus(x[["native"]])$cell_id, ids)
        expect_identical(x[["native"]], p)
    }
    x <- wp03_set()
    x@segmentations$native@universe_fingerprint <- paste(rep("0", 64L), collapse = "")
    expect_error(x[["native"]], "binding")
    x <- wp03_set()
    x@segmentations$native@cell_index[1L] <- 2L
    expect_error(x[["native"]], "fingerprint")
})
