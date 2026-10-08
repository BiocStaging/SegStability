test_that("replicate input selectors fail loudly without convenience coercion", {
    x <- wp10_set()
    for (nm in list(character(), "run_a", c("run_a", "run_a"),
        c("run_a", NA_character_), c("run_a", ""), c("run_a", "absent"),
        1:2, factor(c("run_a", "run_b")), list("run_a", "run_b"),
        matrix(c("run_a", "run_b")), structure(c("run_a", "run_b"), class = "custom"))) {
        expect_error(replicateAgreement(x, nm), "segmentations")
    }
    for (object in list(x@segmentations, x@segmentations[[1L]],
                        compareSegmentations(x, "run_a", "run_b"), NULL)) {
        expect_error(replicateAgreement(object, c("run_a", "run_b")), "SegmentationSet")
    }
    expect_error(replicateAgreement(x), "segmentations")
})

test_that("replicate analysis checks scientific structure universe and assignment integrity", {
    x <- wp10_set()
    bad <- x
    bad@segmentations$run_a@cell_index[1L] <- 2L
    expect_error(replicateAgreement(bad, c("run_a", "run_b")), "fingerprint")
    bad <- x
    bad@molecules$molecule_id[1L] <- "corrupt"
    expect_error(replicateAgreement(bad, c("run_a", "run_b")), "fingerprint")
    bad <- x
    bad@segmentations$run_a@universe_fingerprint <- strrep("b", 64L)
    expect_error(replicateAgreement(bad, c("run_a", "run_b")), "binding")
    bad <- x
    bad@segmentations$run_a@status[1L] <- 4L
    expect_error(replicateAgreement(bad, c("run_a", "run_b")), "status")
    bad <- x
    bad@specimen_id <- NA_character_
    expect_error(replicateAgreement(bad, c("run_a", "run_b")), "specimen_id")
    bad <- x
    names(bad@segmentations) <- c("run_a", "run_a")
    expect_error(replicateAgreement(bad, c("run_a", "run_b")))
    bad <- x
    bad@molecules$feature_index[1L] <- 2L
    expect_error(replicateAgreement(bad, c("run_a", "run_b")), "feature")
})

test_that("integrity is checked once per universe and once per selected run", {
    x <- wp10_set(list(a = c(1L, 1L, 2L, 2L), b = c(1L, 1L, 2L, 2L),
                      c = c(1L, 1L, 2L, 2L), unused = c(1L, 1L, 2L, 2L)))
    counts <- new.env(parent = emptyenv())
    counts$universe <- counts$assignment <- 0L
    trace(".universeFingerprint", where = asNamespace("SegStability"), print = FALSE,
        tracer = function() counts$universe <- counts$universe + 1L)
    trace(".assignmentFingerprint", where = asNamespace("SegStability"), print = FALSE,
        tracer = function() counts$assignment <- counts$assignment + 1L)
    on.exit({
        untrace(".universeFingerprint", where = asNamespace("SegStability"))
        untrace(".assignmentFingerprint", where = asNamespace("SegStability"))
    })
    replicateAgreement(x, c("a", "b", "c"))
    expect_identical(counts$universe, 1L)
    expect_identical(counts$assignment, 3L)
})

test_that("replicate analysis uses existing indices without invoking derived modules", {
    x <- wp10_set()
    stop_dependency <- function(...) stop("Unexpected derived dependency")
    local_mocked_bindings(stateStability = stop_dependency,
        supportSensitivity = stop_dependency, anchorTopology = stop_dependency,
        compareGeometry = stop_dependency, .compareAlignedPartitions = stop_dependency,
        .matchIds = stop_dependency, .alignMolecules = stop_dependency,
        .package = "SegStability")
    expect_identical(replicateAgreement(x, c("run_a", "run_b"))$pairwise_summary$
        ari_jointly_assigned, 1)
})
