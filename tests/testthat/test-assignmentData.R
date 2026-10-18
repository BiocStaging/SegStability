test_that("assignmentData reconstructs canonical typed membership without mutation", {
    for (mkind in c("integer", "character", "LLint", "integer64")) {
        if (mkind == "integer64" && !requireNamespace("bit64", quietly = TRUE)) next
        for (ckind in c("integer", "character", "LLint", "integer64")) {
            if (ckind == "integer64" && !requireNamespace("bit64", quietly = TRUE)) next
            mid <- id_lab_fixture(mkind)$x
            cid <- id_lab_fixture(ckind)$x
            n <- length(mid)
            m <- S4Vectors::DataFrame(molecule_id = mid,
                feature_id = rep("g", n))
            status <- rep(c("assigned", "unassigned", "ambiguous"), length.out = n)
            # Typed missing values without LLint missing subscripts/replacement.
            missing <- switch(ckind, integer = NA_integer_, character = NA_character_,
                LLint = S4Vectors::as.LLint(NA_character_),
                integer64 = bit64::as.integer64(NA_character_))
            pool <- c(cid, missing)
            index <- rep(length(pool), n)
            index[status == "assigned"] <- rep(c(1L, 2L), length.out = sum(status == "assigned"))
            expected <- pool[index]
            order <- rev(seq_len(n))
            p <- SegmentationPartition(m, "native", S4Vectors::DataFrame(cell_id = cid),
                mid[order], status[order], expected[order])
            x <- SegmentationSet(m, list(native = p), "specimen")
            before <- serialize(x, NULL)
            result <- assignmentData(x, "native")
            expect_s4_class(result, "DataFrame")
            expect_identical(colnames(result), c("molecule_id", "status", "cell_id"))
            expect_identical(nrow(result), n)
            expect_id_lab_exact(result$molecule_id, mid)
            expect_identical(result$status, status)
            expect_id_lab_exact(result$cell_id, expected)
            f <- tempfile()
            saveRDS(x, f)
            expect_identical(assignmentData(readRDS(f), "native"), result)
            unlink(f)
            result$status[] <- "unassigned"
            result$molecule_id <- NULL
            expect_identical(serialize(x, NULL), before)
        }
    }
})

test_that("assignmentData preserves ambiguity evidence and zero populations", {
    x <- wp03_set(segmentations = list(native = wp02_partition(ambiguity = wp02_ambiguity())))
    before <- ambiguityData(x[["native"]])
    expect_identical(assignmentData(x, "native")$cell_id, c("c1", NA_character_, NA_character_))
    expect_identical(ambiguityData(x[["native"]]), before)
    for (kind in c("integer", "character", "LLint", "integer64")) {
        if (kind == "integer64" && !requireNamespace("bit64", quietly = TRUE)) next
        ids <- id_lab_fixture(kind)$x[integer()]
        for (n in c(0L, 2L)) {
            m <- S4Vectors::DataFrame(molecule_id = seq_len(n), feature_id = rep("g", n))
            p <- SegmentationPartition(m, "empty", S4Vectors::DataFrame(cell_id = ids),
                seq_len(n), rep("unassigned", n))
            x <- SegmentationSet(m, list(empty = p), "specimen")
            r <- assignmentData(x, "empty")
            expect_identical(nrow(r), n)
            expect_identical(class(r$cell_id), class(ids))
            expect_true(all(is.na(r$cell_id)))
            expect_identical(r$status, rep("unassigned", n))
        }
    }
})

test_that("assignmentData rejects invalid selectors and corrupt binding", {
    x <- wp03_set()
    expect_error(assignmentData(x), "segmentation")
    for (s in list("nat", "absent", "", character(), NA_character_, 1L,
                   c("native", "other"))) expect_error(assignmentData(x, s), "segmentation")
    x@segmentations$native@universe_fingerprint <- paste(rep("0", 64L), collapse = "")
    expect_error(assignmentData(x, "native"), "binding")
    x <- wp03_set()
    x@segmentations$native@cell_index[1L] <- 2L
    expect_error(assignmentData(x, "native"), "fingerprint")
    x <- wp03_set()
    x@molecules$molecule_id[1L] <- "changed"
    expect_error(assignmentData(x, "native"), "fingerprint")
})
