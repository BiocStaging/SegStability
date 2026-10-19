# Assisted-by: OpenAI Codex
test_that("historical contingency fixtures retain graph and dominance semantics", {
    # Validated matrices from 136/144; topology expectations are not ported.
    cases <- list(identity = diag(c(60L, 40L)), split = matrix(c(60L, 40L), 1L),
        merge = matrix(c(60L, 40L), 2L),
        complex = matrix(c(50L, 20L, 30L, 0L), 2L),
        tied = matrix(50L, 2L, 2L), one_tie = matrix(50L, 1L, 2L))
    classes <- list(c("one_to_one", "one_to_one"), "split", "merge",
                    "many_to_many", "many_to_many", "split")
    reciprocal <- c(2L, 1L, 1L, 1L, 0L, 0L)
    for (k in seq_along(cases)) {
        m <- cases[[k]]
        x <- wp04_compare(m)
        h <- correspondenceEdges(x)
        q <- S4Vectors::from(h)
        s <- S4Vectors::to(h)
        e <- S4Vectors::mcols(h)
        expect_identical(e$weight, as.numeric(m[cbind(q, s)]))
        expect_identical(order(q, s, method = "radix"), seq_along(q))
        expect_identical(componentSummary(x)$component_class, classes[[k]])
        expect_equal(sum(e$reciprocal_dominant), reciprocal[k])
        for (side in 1:2) {
            c <- cellSummary(x, c("left", "right")[side])
            index <- if (side == 1L) q else s
            f <- e[[paste0("fraction_", side)]]
            d <- e[[paste0("dominant_", side)]]
            expect_equal(f, e$weight / c$assigned_mass[index])
            for (i in seq_len(nrow(c))) {
                z <- which(index == i)
                expect_identical(d[z], f[z] == max(f[z]))
                expect_identical(d[z], e$weight[z] == max(e$weight[z]))
                expect_equal(c$dominance_count[i], sum(d[z]))
            }
        }
        expect_true(validObject(x))
    }
})

test_that("each assigned/assigned molecule contributes once before aggregation", {
    ss <- wp04_set(matrix(c(3L, 2L, 1L, 0L), 2L), extra = list(
        i = c(1L, 1L, NA_integer_, NA_integer_),
        j = c(NA_integer_, NA_integer_, 2L, NA_integer_),
        s1 = c("assigned", "assigned", "ambiguous", "unassigned"),
        s2 = c("ambiguous", "unassigned", "assigned", "ambiguous")))
    p <- ss@segmentations[[1L]]
    r <- ss@segmentations[[2L]]
    raw <- SegStability:::.sharedPairs(p, r)
    both <- which(p@status == 1L & r@status == 1L)
    expect_identical(raw$molecule_index, both)
    expect_identical(raw$query, p@cell_index[both])
    expect_identical(raw$subject, r@cell_index[both])
    expect_identical(raw$weight, rep(1, length(both)))
    expect_false(anyDuplicated(raw$molecule_index) > 0L)
    x <- compareSegmentations(ss, "left", "right")
    e <- S4Vectors::mcols(correspondenceEdges(x))
    a <- assignmentSummary(x)
    n <- a$n_molecules[a$status_1 == "assigned" & a$status_2 == "assigned"]
    expect_equal(n, length(both))
    expect_lt(nrow(e), nrow(raw))
    expect_equal(sum(raw$weight), sum(e$weight))
    expect_equal(sum(e$weight), n)
    expect_equal(sum(a$n_molecules), x@universe_size)
    for (side in 1:2) {
        c <- cellSummary(x, c("left", "right")[side])
        expect_equal(sum(c$matched_mass), n)
        expect_equal(c$assigned_mass, c$matched_mass + c$unmatched_mass)
        expect_equal(sum(c$assigned_mass), sum(ss@segmentations[[side]]@status == 1L))
        expect_true(any(c$assigned_mass > c$matched_mass))
        expect_equal(c$assigned_mass,
            as.numeric(tabulate(ss@segmentations[[side]]@cell_index, nrow(c))))
        h <- correspondenceEdges(x)
        idx <- if (side == 1L) S4Vectors::from(h) else S4Vectors::to(h)
        for (i in seq_len(nrow(c))) {
            z <- which(idx == i)
            expect_identical(e[[paste0("dominant_", side)]][z],
                e[[paste0("fraction_", side)]][z] ==
                max(e[[paste0("fraction_", side)]][z]))
        }
    }
})

test_that("all status combinations and full census isolates remain explicit", {
    statuses <- c("assigned", "unassigned", "ambiguous")
    extra <- list(s1 = rep(statuses, each = 3L), s2 = rep(statuses, 3L))
    extra$i <- ifelse(extra$s1 == "assigned", 1L, NA_integer_)
    extra$j <- ifelse(extra$s2 == "assigned", 1L, NA_integer_)
    x <- wp04_compare(matrix(0L, 2L, 2L), extra = extra)
    expect_identical(assignmentSummary(x)$n_molecules, rep(1, 9L))
    expect_identical(cellSummary(x, "left")$assigned_mass, c(3, 0))
    expect_identical(cellSummary(x, "left")$matched_fraction, c(1/3, NA_real_))
    expect_identical(componentSummary(x)$component_id, 1:3)
    expect_identical(cellSummary(x, "right")$component_id, c(1L, 3L))
    expect_identical(componentSummary(x)$component_class,
                     c("one_to_one", "unmatched", "unmatched"))
    z <- wp04_compare(matrix(0L, 2L, 2L), extra = list(i = 1L,
        j = NA_integer_, s1 = "assigned", s2 = "unassigned"))
    expect_length(correspondenceEdges(z), 0L)
    expect_identical(cellSummary(z, "left")$assigned_mass, c(1, 0))
    expect_identical(cellSummary(z, "left")$dominance_count, c(0L, 0L))
    expect_identical(cellSummary(z, "left")$matched_fraction, c(0, NA_real_))
    expect_identical(componentSummary(z)$component_class, rep("unmatched", 4L))
    for (dims in list(c(0L, 0L), c(0L, 2L), c(2L, 0L))) {
        y <- wp04_compare(matrix(0L, dims[1L], dims[2L]))
        expect_true(validObject(y))
        expect_equal(nrow(componentSummary(y)), sum(dims))
        expect_type(componentSummary(y)$component_class, "character")
        expect_equal(S4Vectors::nLnode(correspondenceEdges(y)), dims[1L])
        expect_equal(S4Vectors::nRnode(correspondenceEdges(y)), dims[2L])
    }
})

test_that("a weak positive edge joins components before any selection", {
    x <- wp04_compare(matrix(c(9L, 0L, 1L, 8L), 2L))
    h <- correspondenceEdges(x)
    expect_length(h, 3L)
    expect_identical(componentSummary(x)$component_class, "many_to_many")
    expect_identical(cellSummary(x, "left")$component_id, c(1L, 1L))
    expect_identical(S4Vectors::mcols(h)$dominant_1, c(TRUE, FALSE, TRUE))
    expect_identical(S4Vectors::mcols(h)$reciprocal_dominant,
                     c(TRUE, FALSE, TRUE))
    rev <- compareSegmentations(wp04_set(matrix(c(9L, 0L, 1L, 8L), 2L)),
                                "right", "left")
    expect_equal(sum(S4Vectors::mcols(correspondenceEdges(rev))$weight), 18)
    expect_identical(cellSummary(x, "left")$assigned_mass,
                     cellSummary(rev, "left")$assigned_mass)
})

test_that("ties remain fraction maxima with mass unresolved on the other side", {
    x <- wp04_compare(matrix(2L, 1L, 2L), extra = list(
        i = c(1L, 1L, NA_integer_), j = c(NA_integer_, NA_integer_, 1L),
        s1 = c("assigned", "assigned", "unassigned"),
        s2 = c("ambiguous", "unassigned", "assigned")))
    e <- S4Vectors::mcols(correspondenceEdges(x))
    expect_identical(e$weight, c(2, 2))
    expect_identical(e$fraction_1, c(1/3, 1/3))
    expect_identical(e$fraction_2, c(2/3, 1))
    expect_identical(e$dominant_1, c(TRUE, TRUE))
    expect_identical(e$unique_dominant_1, c(FALSE, FALSE))
    expect_identical(e$unique_dominant_2, c(TRUE, TRUE))
    expect_identical(e$reciprocal_dominant, c(FALSE, FALSE))
    expect_identical(cellSummary(x, "left")$dominance_count, 2L)
    expect_identical(e$dominant_1, e$fraction_1 == max(e$fraction_1))
})

test_that("historical empty opposite census preserves positive own denominators", {
    x <- wp04_compare(matrix(0L, 2L, 0L), extra = list(
        i = 1:2, j = c(NA_integer_, NA_integer_),
        s1 = c("assigned", "assigned"), s2 = c("unassigned", "ambiguous")))
    expect_identical(cellSummary(x, "left")$assigned_mass, c(1, 1))
    expect_identical(cellSummary(x, "left")$matched_mass, c(0, 0))
    expect_identical(cellSummary(x, "left")$unmatched_mass, c(1, 1))
    expect_identical(cellSummary(x, "left")$matched_fraction, c(0, 0))
    expect_identical(nrow(cellSummary(x, "right")), 0L)
    expect_identical(componentSummary(x)$component_class,
                     c("unmatched", "unmatched"))
    expect_identical(sum(assignmentSummary(x)$n_molecules), 2)
})

test_that("comparison preserves exact identifier families and source metadata", {
    trace(".from_LLint_to_numeric", where = asNamespace("S4Vectors"),
          tracer = quote(stop("unsafe double conversion")), print = FALSE)
    on.exit(untrace(".from_LLint_to_numeric",
                    where = asNamespace("S4Vectors")), add = TRUE)
    if (requireNamespace("bit64", quietly = TRUE)) {
        trace("as.double.integer64", where = asNamespace("bit64"),
              tracer = quote(stop("unsafe double conversion")), print = FALSE)
        on.exit(untrace("as.double.integer64", where = asNamespace("bit64")),
                add = TRUE)
    }
    for (kind in c("integer", "character", "LLint", "integer64")) {
        if (kind == "integer64" && !requireNamespace("bit64", quietly = TRUE)) {
            next
        }
        ids <- id_lab_fixture(kind)$x
        m <- S4Vectors::DataFrame(molecule_id = ids,
                                 feature_id = rep("g", length(ids)))
        run <- S4Vectors::SimpleList(nested = list(seed = 7L),
            missing = NA_character_, raw = as.raw(c(0, 255)))
        p <- SegmentationPartition(m, "left",
            S4Vectors::DataFrame(cell_id = ids), ids,
            rep("assigned", length(ids)), ids, run_metadata = run)
        # Reversed source molecule order must align before comparison.
        q <- SegmentationPartition(m, "right",
            S4Vectors::DataFrame(cell_id = ids), rev(ids),
            rep("assigned", length(ids)), ids)
        ss <- SegmentationSet(m, list(left = p, right = q), "specimen",
            coordinate_info = S4Vectors::SimpleList(units = "unknown"))
        x <- compareSegmentations(ss, "left", "right")
        expect_id_lab_exact(cellSummary(x, "left")$cell_id, ids)
        expect_id_lab_exact(cellSummary(x, "right")$cell_id, ids)
        expect_identical(S4Vectors::to(correspondenceEdges(x)),
                         rev(seq_along(ids)))
        expect_identical(x@coordinate_info, ss@coordinate_info)
        expect_identical(x@metadata$run_metadata$left, run)
        expect_true(is.na(provenance(x)$source))
        file <- tempfile(fileext = ".rds")
        saveRDS(x, file)
        y <- readRDS(file)
        unlink(file)
        expect_identical(y, x)
        expect_true(validObject(y))
        expect_identical(compareSegmentations(ss, "left", "right"), x)
    }
})

test_that("comparison hashes the universe once and does not hash in validity", {
    ss <- wp04_set()
    count <- new.env(parent = emptyenv())
    count$n <- 0L
    trace(".universeFingerprint", where = asNamespace("SegStability"),
        print = FALSE, tracer = function() count$n <- count$n + 1L)
    on.exit(untrace(".universeFingerprint", where = asNamespace("SegStability")))
    x <- compareSegmentations(ss, "left", "right")
    expect_identical(count$n, 1L)
    expect_true(validObject(x))
    expect_identical(count$n, 1L)
    bad <- ss
    bad@molecules$molecule_id[1L] <- "corrupt"
    expect_true(validObject(bad))
    expect_error(compareSegmentations(bad, "left", "right"),
                 "Universe fingerprint")
})

test_that("labels and census permutations do not establish correspondence", {
    m <- matrix(c(3L, 1L, 2L, 0L), 2L)
    x <- wp04_compare(m, ids1 = c("same", "other"),
                       ids2 = c("other", "same"))
    y <- wp04_compare(m, ids1 = c("renamed1", "renamed2"),
                       ids2 = c("different1", "different2"))
    expect_identical(correspondenceEdges(x), correspondenceEdges(y))
    z <- wp04_compare(m[2:1, 2:1], ids1 = c("other", "same"),
                       ids2 = c("same", "other"))
    keyed <- function(cmp) {
        h <- correspondenceEdges(cmp)
        d <- data.frame(left = cellSummary(cmp, "left")$cell_id[S4Vectors::from(h)],
            right = cellSummary(cmp, "right")$cell_id[S4Vectors::to(h)],
            weight = S4Vectors::mcols(h)$weight)
        d <- d[order(d$left, d$right, method = "radix"), ]
        rownames(d) <- NULL
        d
    }
    expect_identical(keyed(x), keyed(z))
    expect_equal(S4Vectors::mcols(correspondenceEdges(x))$weight, c(3, 2, 1))
})

test_that("comparison validates API, bindings and source integrity", {
    ss <- wp04_set()
    expect_identical(names(formals(compareSegmentations)),
                     c("object", "segmentation1", "segmentation2"))
    expect_error(compareSegmentations(ss, "left", "left"), "distinct")
    expect_error(compareSegmentations(ss, "left", "missing"), "segmentation")
    expect_error(compareSegmentations(ss, NA_character_, "right"), "segmentation")
    expect_error(compareSegmentations(ss, 1L, "right"), "segmentation")
    expect_error(compareSegmentations(ss@segmentations[[1L]], "left", "right"),
                 "SegmentationSet")
    bad <- ss
    bad@segmentations[[1L]]@cell_index[1L] <- 1L
    expect_error(compareSegmentations(bad, "left", "right"), "fingerprint")
    bad <- ss
    bad@molecules$feature_index[1L] <- 2L
    expect_error(compareSegmentations(bad, "left", "right"))
    bad <- ss
    bad@segmentations[[1L]]@universe_fingerprint <- strrep("b", 64L)
    expect_error(compareSegmentations(bad, "left", "right"), "binding")
    x <- compareSegmentations(ss, "left", "right")
    expect_identical(x@assignment_fingerprints,
        vapply(ss@segmentations, function(p) p@assignment_fingerprint, character(1)))
    expect_identical(provenance(x)$set, provenance(ss))
    expect_identical(provenance(x)$segmentations$left,
                     provenance(ss@segmentations[[1L]]))
    expect_identical(ss, wp04_set())
})
