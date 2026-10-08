test_that("replicateAgreement exposes exactly the approved signature and schema", {
    expect_true("replicateAgreement" %in% getNamespaceExports("SegStability"))
    expect_identical(names(formals(replicateAgreement)), c("object", "segmentations"))
    x <- wp10_set()
    r <- replicateAgreement(x, c("run_a", "run_b"))
    expect_s4_class(r, "SimpleList")
    expect_identical(names(r), c("pairwise_summary", "assignment_status_contingency"))
    expect_s4_class(r$pairwise_summary, "DataFrame")
    expect_identical(colnames(r$pairwise_summary), c("segmentation_1", "segmentation_2",
        "n_molecules", "n_status_agree", "n_status_disagree", "status_agreement_fraction",
        "n_jointly_assigned", "jointly_assigned_fraction", "ari_jointly_assigned"))
    expect_s4_class(r$assignment_status_contingency, "DataFrame")
    expect_identical(colnames(r$assignment_status_contingency),
        c("segmentation_1", "segmentation_2", "status_1", "status_2", "n_molecules"))
    expect_identical(r$pairwise_summary$ari_jointly_assigned, 1)
    expect_identical(r$pairwise_summary$n_jointly_assigned, 4)
    expect_identical(S4Vectors::metadata(r)$segmentations, c("run_a", "run_b"))
})

test_that("all status combinations are explicit and conserve the full universe", {
    left <- c(1L, 1L, 1L, 0L, 0L, 0L, -1L, -1L, -1L)
    right <- rep(c(1L, 0L, -1L), 3L)
    r <- replicateAgreement(wp10_set(list(a = left, b = right)), c("a", "b"))
    status <- c("assigned", "unassigned", "ambiguous")
    tab <- r$assignment_status_contingency
    expect_identical(tab$status_1, rep(status, each = 3L))
    expect_identical(tab$status_2, rep(status, 3L))
    expect_identical(tab$n_molecules, rep(1, 9L))
    s <- r$pairwise_summary
    expect_identical(s$n_molecules, 9)
    expect_identical(s$n_status_agree, 3)
    expect_identical(s$n_status_disagree, 6)
    expect_identical(s$status_agreement_fraction, 1 / 3)
    expect_identical(s$n_jointly_assigned, 1)
    expect_identical(s$jointly_assigned_fraction, 1 / 9)
    expect_identical(s$ari_jointly_assigned, NA_real_)
})

test_that("unresolved populations and degenerate partitions never acquire an ARI", {
    for (a in list(integer(), rep(0L, 4L), rep(-1L, 4L), 1L,
                   rep(1L, 4L), 1:4)) {
        x <- wp10_set(list(a = a, b = a))
        before <- serialize(x, NULL)
        r <- replicateAgreement(x, c("a", "b"))
        s <- r$pairwise_summary
        expect_identical(s$ari_jointly_assigned, NA_real_)
        expect_identical(sum(r$assignment_status_contingency$n_molecules),
                         as.numeric(length(a)))
        expect_identical(s$n_status_disagree, 0)
        expect_identical(s$status_agreement_fraction, if (length(a)) 1 else NA_real_)
        expect_identical(s$jointly_assigned_fraction,
            if (length(a)) sum(a > 0L) / length(a) else NA_real_)
        expect_identical(serialize(x, NULL), before)
    }
    r <- replicateAgreement(wp10_set(list(a = rep(0L, 4L), b = rep(-1L, 4L))),
                             c("a", "b"))
    expect_identical(r$pairwise_summary$n_status_disagree, 4)
    expect_identical(r$pairwise_summary$n_jointly_assigned, 0)
})

test_that("partition changes match an independent molecule-pair oracle", {
    cases <- list(
        identity = list(c(1L, 1L, 2L, 2L), c(2L, 2L, 1L, 1L)),
        reassignment = list(c(1L, 1L, 2L, 2L), c(1L, 2L, 2L, 2L)),
        split = list(c(1L, 1L, 2L, 2L), c(1L, 2L, 3L, 3L)),
        merge = list(c(1L, 2L, 3L, 3L), c(1L, 1L, 2L, 2L)),
        complex_ties = list(c(1L, 1L, 2L, 2L), c(1L, 2L, 1L, 2L)),
        zero_correspondence = list(c(1L, 1L, 2L, 2L, 3L, 3L),
                                   c(1L, 1L, 2L, 2L, 0L, 0L)),
        unresolved = list(c(1L, 1L, 2L, 2L, 1L, 0L, -1L),
                          c(1L, 1L, 2L, 2L, 0L, -1L, 1L)))
    for (case in cases) {
        names(case) <- c("a", "b")
        x <- wp10_set(case)
        s <- replicateAgreement(x, c("a", "b"))$pairwise_summary
        expect_equal(s$ari_jointly_assigned, wp10_pair_oracle(case$a, case$b))
        expect_identical(s$n_jointly_assigned, as.numeric(sum(case$a > 0 & case$b > 0)))
        expect_identical(s$n_status_agree + s$n_status_disagree, s$n_molecules)
        expect_identical(nrow(cellCensus(x@segmentations$a)), max(case$a) + 1L)
    }
})

test_that("multi-run pair order follows the supplied selectors and is symmetric", {
    x <- wp10_set(list(a = c(1L, 1L, 2L, 2L), b = c(1L, 2L, 1L, 2L),
                      c = c(1L, 1L, 1L, 2L), d = c(1L, 0L, 2L, -1L)))
    r <- replicateAgreement(x, c("d", "b", "a", "c"))
    s <- r$pairwise_summary
    expect_identical(s$segmentation_1, c("d", "d", "d", "b", "b", "a"))
    expect_identical(s$segmentation_2, c("b", "a", "c", "a", "c", "c"))
    expect_identical(nrow(r$assignment_status_contingency), 54L)
    for (i in seq_len(nrow(s))) {
        nm <- c(s$segmentation_1[i], s$segmentation_2[i])
        f <- replicateAgreement(x, nm)
        b <- replicateAgreement(x, rev(nm))
        expect_identical(f$pairwise_summary[, -(1:2)], b$pairwise_summary[, -(1:2)])
        expect_identical(f$assignment_status_contingency$n_molecules,
            b$assignment_status_contingency$n_molecules[c(1, 4, 7, 2, 5, 8, 3, 6, 9)])
    }
    expect_identical(r, replicateAgreement(x, c("d", "b", "a", "c")))
})

test_that("run metadata inputs and serialized results remain unchanged", {
    metadata <- list(S4Vectors::SimpleList(method = "unknown", seed = NA_integer_,
        run_id = "local", nested = S4Vectors::SimpleList(bytes = as.raw(c(0, 255)))),
        S4Vectors::SimpleList(method = "method", seed = 7L, run_id = "local",
            table = S4Vectors::DataFrame(value = c("x", NA_character_))))
    x <- wp10_set(run_metadata = metadata)
    original <- serialize(x, NULL)
    r <- replicateAgreement(x, c("run_b", "run_a"))
    expect_identical(serialize(x, NULL), original)
    expect_identical(names(S4Vectors::metadata(r)), c("segmentations", "run_metadata"))
    expect_identical(S4Vectors::metadata(r)$run_metadata$run_a, metadata[[1L]])
    expect_identical(S4Vectors::metadata(r)$run_metadata$run_b, metadata[[2L]])
    path <- tempfile()
    on.exit(unlink(path))
    saveRDS(list(x = x, r = r), path)
    restored <- readRDS(path)
    expect_identical(restored$r, r)
    expect_identical(replicateAgreement(restored$x, c("run_b", "run_a")), r)
})
