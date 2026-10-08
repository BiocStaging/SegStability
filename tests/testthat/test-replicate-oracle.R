test_that("sparse ARI matches independent tiny pair enumeration", {
    configurations <- as.matrix(expand.grid(rep(list(1:2), 4L)))
    for (i in seq_len(nrow(configurations))) {
        left <- unname(configurations[i, ])
        for (j in seq_len(nrow(configurations))) {
            right <- unname(configurations[j, ])
            pairs <- S4Vectors::DataFrame(query = left, subject = right)
            h <- .aggregatePairs(pairs, 2L, 2L)
            expect_equal(.replicateARI(h), wp10_pair_oracle(left, right), tolerance = 1e-14)
        }
    }
})

test_that("validated historical sparse contingency arithmetic reproduces saved ARI", {
    # Independent exact-integer aggregates of the three read-only ROI1 tables.
    # The producing script and independent receipt are recorded in the local audit.
    values <- list(
        c("27107867", "28189168", "28245360", "19163631651", "0.9606257741452051"),
        c("27202443", "28237319", "28305123", "19120781235", "0.9621395697764434"),
        c("27102690", "28201917", "28200823", "19143472285", "0.9609841324123201"))
    for (v in values) {
        expect_equal(.replicateARIFromCounts(wp10_ll(v[1L]), wp10_ll(v[2L]),
            wp10_ll(v[3L]), wp10_ll(v[4L])), as.numeric(v[5L]), tolerance = 1e-14)
    }
})
