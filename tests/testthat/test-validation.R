# Assisted-by: OpenAI Codex
test_that("reports have the exact schema and strict error policy", {
    p <- wp02_partition()
    before <- serialize(p, NULL)
    r <- validateSegmentation(p, strict = TRUE)
    expect_s4_class(r, "DataFrame")
    expect_identical(colnames(r), c("check", "status", "severity", "message"))
    expect_true(all(vapply(as.list(r), is.character, logical(1))))
    expect_false(validation_defect(r))
    expect_identical(validation_row(r, "integrity.assignment")$status, "PASS")
    expect_identical(validation_row(r, "universe.verification")$status, "UNKNOWN")
    expect_identical(validation_row(validateSegmentation(p, check_integrity = FALSE),
                                   "integrity.assignment")$status, "UNKNOWN")
    expect_identical(r, validateSegmentation(p))
    expect_identical(serialize(p, NULL), before)
    expect_identical(unserialize(serialize(r, NULL)), r)
    p@status[1L] <- 9L
    expect_true(validation_defect(validateSegmentation(p)))
    expect_error(validateSegmentation(p, strict = TRUE), "Validation failed")
    expect_true(validation_defect(validateSegmentation(list())))
})

test_that("partition diagnostics catch corruption and preserve scientific states", {
    p <- wp02_partition(ambiguity = wp02_ambiguity())
    changes <- list(function(x) { x@cells$cell_id[2L] <- "c1"; x },
        function(x) { x@cell_index[3L] <- 1L; x },
        function(x) { x@ambiguity$evidence$candidate_count <- 8L; x },
        function(x) { x@ambiguity$candidates <- S4Vectors::Hits(3L, 9L, 3L, 9L); x },
        function(x) { x@assignment_fingerprint <- strrep("a", 64); x })
    for (change in changes) {
        bad <- change(p)
        before <- serialize(bad, NULL)
        expect_true(validation_defect(validateSegmentation(bad)))
        expect_identical(serialize(bad, NULL), before)
    }
    expect_false(validation_defect(validateSegmentation(p)))
    expect_identical(validation_row(validateSegmentation(wp02_partition()),
                                   "ambiguity.completeness")$status, "UNKNOWN")
    m <- wp02_molecules()[integer(), ]
    empty <- SegmentationPartition(m, "empty",
        S4Vectors::DataFrame(cell_id = character()), character(), character())
    expect_false(validation_defect(validateSegmentation(empty)))
})

test_that("set integrity verifies universe once and binds every partition", {
    x <- wp03_set()
    before <- serialize(x, NULL)
    r <- validateSegmentationSet(x, strict = TRUE)
    expect_false(validation_defect(r))
    expect_identical(validation_row(r, "integrity.universe")$status, "PASS")
    expect_identical(serialize(x, NULL), before)
    bad <- x
    bad@molecules$molecule_id[1L] <- "same-size-wrong-universe"
    expect_true(validation_defect(validateSegmentationSet(bad)))
    expect_false(validation_defect(validateSegmentationSet(bad,
                                                          check_integrity = FALSE)))
    for (change in list(function(x) { names(x@segmentations) <- c("a", "a"); x },
        function(x) { x@molecules$feature_index[1L] <- 99L; x },
        function(x) { x@specimen_id <- NA_character_; x },
        function(x) { x@segmentations[[1L]]@status <- 9L; x })) {
        expect_true(validation_defect(validateSegmentationSet(change(x))))
    }
})

test_that("coordinate measurements accept finite doubles and align only by ID", {
    d <- S4Vectors::DataFrame(molecule_id = c("m1", "m2", "m3"),
        x = c(0.25, 2.5, -1.5), y = c(2, 4, 6))
    before <- serialize(d, NULL)
    expect_false(validation_defect(validateCoordinates(d, strict = TRUE)))
    expect_identical(validation_row(validateCoordinates(NULL),
                                   "coordinates.presence")$status, "NOT_APPLICABLE")
    for (v in list(NA_real_, NaN, Inf, "bad", TRUE)) {
        bad <- d
        bad$x <- rep(v, 3)
        expect_true(validation_defect(validateCoordinates(bad)))
    }
    expect_true(validation_defect(validateCoordinates(d[, c("molecule_id", "x")])))
    bad <- d
    bad$molecule_id[2L] <- "m1"
    expect_true(validation_defect(validateCoordinates(bad)))
    m <- wp02_molecules()
    m$x <- d$x
    m$y <- d$y
    p <- SegmentationPartition(m, "s", S4Vectors::DataFrame(cell_id = "a"),
        m$molecule_id, rep("assigned", 3), rep("a", 3))
    ref <- SegmentationSet(m, list(s = p), "specimen")
    expect_false(validation_defect(validateCoordinates(d[c(3L, 1L, 2L), ],
                                                       reference = ref)))
    bad <- d
    bad$x[1L] <- 9
    expect_true(validation_defect(validateCoordinates(bad, reference = ref)))
    expect_identical(serialize(d, NULL), before)
})

test_that("provenance completeness remains distinct from authenticity", {
    p <- wp02_partition()
    r <- validateProvenance(p, strict = TRUE)
    expect_false(validation_defect(r))
    expect_true(any(r$status == "UNKNOWN"))
    expect_identical(validation_row(r, "provenance.authenticity")$status, "UNKNOWN")
    expect_true(validation_defect(validateProvenance(S4Vectors::SimpleList())))
    pr <- provenance(p)
    pr$threads <- "eight"
    expect_true(validation_defect(validateProvenance(pr)))
    pr <- provenance(p)
    pr$runtime <- -1
    expect_error(validateProvenance(pr, strict = TRUE), "Validation failed")
})

test_that("spatial declarations diagnose conflicts without certifying registration", {
    d <- S4Vectors::DataFrame(molecule_id = "m", x = 1.5, y = 2.5)
    ci <- S4Vectors::SimpleList(units = "micrometres", origin = c(0, 0),
        axis_directions = c("right", "up"), transforms = S4Vectors::SimpleList())
    r <- validateCoordinates(d, ci)
    expect_false(validation_defect(r))
    expect_identical(validation_row(r, "coordinates.registration_evidence")$status,
                     "UNKNOWN")
})

test_that("known provenance fields and partial evidence have distinct outcomes", {
    pr <- provenance(wp02_partition())
    for (field in names(pr)) {
        pr[[field]] <- switch(field, threads = 2L, runtime = 0.5,
            parameters = S4Vectors::SimpleList(seed = 1L),
            observation_window = c(xmin = 0, xmax = 2, ymin = 0, ymax = 2),
            "known")
    }
    before <- serialize(pr, NULL)
    r <- validateProvenance(pr, strict = TRUE)
    expect_false(validation_defect(r))
    expect_true(all(r$status[r$check != "provenance.authenticity"] == "PASS"))
    expect_identical(serialize(pr, NULL), before)
    for (field in names(pr)) {
        absent <- pr[names(pr) != field]
        row <- validation_row(validateProvenance(absent), paste0("provenance.", field))
        expect_identical(row$status, "FAIL")
        expect_identical(row$severity, "error")
    }
    for (entry in list(c("threads", "numeric"), c("runtime", "infinite"),
                       c("source", "list"), c("observation_window", "bad"))) {
        bad <- pr
        bad[[entry[1L]]] <- switch(entry[2L], numeric = 2, infinite = Inf,
            list = list("source"), bad = c(xmin = 0, xmax = 0, ymin = 0, ymax = 1))
        expect_true(validation_defect(validateProvenance(bad)))
    }
    pr$source_hashes <- S4Vectors::SimpleList(input = NA_character_)
    expect_identical(validation_row(validateProvenance(pr),
                                   "provenance.source_hashes")$status, "UNKNOWN")
    expect_true(validation_defect(validateProvenance(list(source = "x"))))
    for (x in list(wp03_set(), compareSegmentations(wp03_set(), "native", "other"))) {
        expect_false(validation_defect(validateProvenance(x)))
    }
})

test_that("coordinate checks preserve every scientific ID family", {
    ids <- list(c(1L, 2L), c("9223372036854775809", "9223372036854775810"),
        S4Vectors::as.LLint(c("9007199254740992", "9007199254740993")))
    if (requireNamespace("bit64", quietly = TRUE)) {
        ids[[length(ids) + 1L]] <- bit64::as.integer64(c("9007199254740992",
                                                       "9007199254740993"))
    }
    for (id in ids) {
        d <- S4Vectors::DataFrame(molecule_id = id, x = c(1.5, 2.5), y = c(2, 3))
        before <- serialize(d, NULL)
        expect_false(validation_defect(validateCoordinates(d, strict = TRUE)))
        expect_identical(serialize(d, NULL), before)
    }
    d <- S4Vectors::DataFrame(molecule_id = c(1, 2), x = c(1.5, 2.5), y = c(2, 3))
    expect_true(validation_defect(validateCoordinates(d)))
    d$molecule_id <- c("m1", "m2")
    d$x <- matrix(c(1, 2), ncol = 1)
    expect_true(validation_defect(validateCoordinates(d)))
    d$x <- c(1, 2)
    d$z <- c(1, Inf)
    expect_true(validation_defect(validateCoordinates(d)))
    d$z <- c(1, 2)
    d$molecule_id <- c("m1", NA_character_)
    expect_true(validation_defect(validateCoordinates(d)))
    expect_true(validation_defect(validateCoordinates(matrix(1:4, ncol = 2))))
})

test_that("coordinate reference conflicts never infer registration or repair", {
    m <- wp02_molecules()
    m$x <- c(1.5, 2.5, 3.5)
    m$y <- c(2, 3, 4)
    p <- SegmentationPartition(m, "s", S4Vectors::DataFrame(cell_id = "c"),
        m$molecule_id, rep("assigned", 3), rep("c", 3))
    ci <- S4Vectors::SimpleList(units = "micrometres", origin = c(0, 0),
        axis_directions = c("right", "up"), transforms = S4Vectors::SimpleList(),
        specimen_id = "specimen", crs = NA_character_, registration = NA_character_)
    x <- SegmentationSet(m, list(s = p), "specimen", coordinate_info = ci)
    before <- serialize(x, NULL)
    expect_false(validation_defect(validateCoordinates(x, reference = x)))
    for (field in c("units", "origin", "axis_directions", "specimen_id")) {
        bad <- ci
        bad[[field]] <- "different"
        expect_true(validation_defect(validateCoordinates(x, bad, reference = x)))
        expect_error(validateCoordinates(x, bad, reference = x, strict = TRUE),
                     "Validation failed")
    }
    d <- moleculeData(x)
    expect_true(validation_defect(validateCoordinates(d[-1L, ], reference = x)))
    extra <- d
    extra$molecule_id[1L] <- "extra"
    expect_true(validation_defect(validateCoordinates(extra, reference = x)))
    expect_true(validation_defect(validateCoordinates(d, reference = list())))
    expect_identical(serialize(x, NULL), before)
})

test_that("nested malformed objects and incomplete evidence stay diagnosable", {
    for (count in c(NA_integer_, 2L)) {
        evidence <- wp02_ambiguity()
        evidence$evidence$candidate_count <- count
        evidence$evidence$candidate_ids_known <- FALSE
        evidence$candidates <- evidence$candidates[integer(), ]
        p <- wp02_partition(ambiguity = evidence)
        r <- validateSegmentation(p)
        expect_false(validation_defect(r))
        expect_identical(validation_row(r, "ambiguity.completeness")$status, "UNKNOWN")
    }
    p@ambiguity$evidence$candidate_ids_known <- TRUE
    expect_true(validation_defect(validateSegmentation(p)))
    x <- wp03_set()
    x@segmentations[[1L]]@universe_fingerprint <- strrep("f", 64)
    expect_true(validation_defect(validateSegmentationSet(x)))
    x <- wp03_set()
    x@segmentations[[1L]] <- list()
    expect_true(validation_defect(validateSegmentationSet(x)))
    x <- wp03_set()
    x@coordinate_info <- S4Vectors::SimpleList(list(a = 1, b = 2))
    names(x@coordinate_info) <- c("a", "a")
    expect_true(validation_defect(validateSegmentationSet(x)))
    expect_error(validateSegmentation(wp02_partition(), check_integrity = NA),
                 "check_integrity")
    expect_error(validateSegmentation(wp02_partition(), strict = NA), "strict")
})

test_that("arbitrary ancillary metadata retains the prior accepted contract", {
    p <- wp02_partition(run_metadata = S4Vectors::SimpleList(7L),
                        metadata = S4Vectors::SimpleList("opaque"))
    x <- wp03_set(metadata = S4Vectors::SimpleList("opaque"))
    expect_false(validation_defect(validateSegmentation(p)))
    expect_false(validation_defect(validateSegmentationSet(x)))
})

test_that("set validation hashes the universe once and assignments independently", {
    x <- wp03_set()
    tally <- new.env(parent = emptyenv())
    tally$universe <- tally$assignment <- 0L
    universe <- .universeFingerprint
    assignment <- .assignmentFingerprint
    local_mocked_bindings(.universeFingerprint = function(...) {
        tally$universe <- tally$universe + 1L
        universe(...)
    }, .assignmentFingerprint = function(...) {
        tally$assignment <- tally$assignment + 1L
        assignment(...)
    })
    expect_false(validation_defect(validateSegmentationSet(x)))
    expect_identical(tally$universe, 1L)
    expect_identical(tally$assignment, 2L)
    validateSegmentationSet(x, check_integrity = FALSE)
    expect_identical(tally$universe, 1L)
    expect_identical(tally$assignment, 2L)
})

test_that("empty populations remain valid set and coordinate inputs", {
    m <- S4Vectors::DataFrame(molecule_id = character(), feature_id = character(),
                              x = numeric(), y = numeric())
    p <- SegmentationPartition(m, "empty", S4Vectors::DataFrame(cell_id = character()),
                               character(), character())
    x <- SegmentationSet(m, list(empty = p), "specimen")
    expect_false(validation_defect(validateSegmentationSet(x, strict = TRUE)))
    expect_false(validation_defect(validateCoordinates(x, reference = x, strict = TRUE)))
})

test_that("census completeness cannot be inferred from observed assignment mass", {
    r <- validateSegmentation(wp02_partition())
    expect_identical(validation_row(r, "census.source_completeness")$status, "UNKNOWN")
})

test_that("core validation retains exact molecule and census ID representations", {
    for (kind in c("integer", "character", "LLint", "integer64")) {
        if (kind == "integer64" && !requireNamespace("bit64", quietly = TRUE)) {
            next
        }
        id <- id_lab_fixture(kind)$x
        m <- S4Vectors::DataFrame(molecule_id = id,
                                  feature_id = rep("g", length(id)))
        p <- SegmentationPartition(m, "s", S4Vectors::DataFrame(cell_id = id),
                                   id, rep("assigned", length(id)), id)
        x <- SegmentationSet(m, list(s = p), "specimen")
        before <- serialize(x, NULL)
        expect_false(validation_defect(validateSegmentation(p, strict = TRUE)))
        expect_false(validation_defect(validateSegmentationSet(x, strict = TRUE)))
        expect_identical(serialize(x, NULL), before)
    }
})

test_that("backend warnings remain inspectable instead of escaping the report", {
    expect_silent(r <- .validationCheck("probe", { warning("backend evidence"); TRUE }))
    expect_identical(r$status, "UNKNOWN")
    expect_identical(r$severity, "warning")
    expect_match(r$message, "backend evidence")
    expect_silent(r <- .validationCheck("probe", {
        warning("backend evidence")
        stop("failed check")
    }, severity = "warning"))
    expect_identical(r$status, "FAIL")
    expect_identical(r$severity, "warning")
    expect_match(r$message, "failed check")
    expect_match(r$message, "backend evidence")
    expect_silent(.validationFinish(r, TRUE))
})

test_that("explicit metadata cannot hide malformed stored or reference declarations", {
    x <- wp03_set()
    x@molecules$x <- c(1.5, 2.5, 3.5)
    x@molecules$y <- c(2.5, 3.5, 4.5)
    x@coordinate_info <- S4Vectors::SimpleList(a = 1, b = 2)
    names(x@coordinate_info) <- c("a", "a")
    r <- validateCoordinates(x, coordinate_info = S4Vectors::SimpleList())
    expect_identical(validation_row(r, "coordinates.declarations.stored")$status, "FAIL")
    r <- validateCoordinates(moleculeData(x), reference = x)
    expect_identical(validation_row(r, "coordinates.reference")$status, "FAIL")
    x@coordinate_info <- S4Vectors::SimpleList(specimen_id = "wrong-specimen")
    r <- validateCoordinates(moleculeData(x), reference = x)
    expect_identical(validation_row(r, "coordinates.reference")$status, "FAIL")
})

test_that("measurement equality does not impose identifier representation rules", {
    m <- wp02_molecules()
    m$x <- 1:3
    m$y <- c(2, 3, 4)
    p <- SegmentationPartition(m, "s", S4Vectors::DataFrame(cell_id = "a"),
        m$molecule_id, rep("assigned", 3), rep("a", 3))
    x <- SegmentationSet(m, list(s = p), "specimen")
    d <- moleculeData(x)
    d$x <- as.double(d$x)
    r <- validateCoordinates(d, reference = x)
    expect_identical(validation_row(r, "coordinates.alignment")$status, "PASS")
})

test_that("a reference specimen never supplies missing external declarations", {
    m <- wp02_molecules()
    m$x <- c(1.5, 2.5, 3.5)
    m$y <- c(2, 3, 4)
    p <- SegmentationPartition(m, "s", S4Vectors::DataFrame(cell_id = "a"),
        m$molecule_id, rep("assigned", 3), rep("a", 3))
    x <- SegmentationSet(m, list(s = p), "specimen")
    r <- validateCoordinates(moleculeData(x), reference = x)
    expect_identical(validation_row(r, "coordinates.specimen_id")$status, "UNKNOWN")
    r <- validateCoordinates(x, reference = x)
    expect_identical(validation_row(r, "coordinates.specimen_id")$status, "PASS")
})
