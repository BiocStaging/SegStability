test_that("matching maximizes actual IoU before applying exact ties", {
    cases <- list(
        c(.5, .5, .5, .5),
        c(.5, .5, .5, .5 + 2^-52),
        c(.9, .8, .8, .1),
        c(1, .5, .5, 0),
        c(.7, .4, .3, .2)
    )
    for (w in cases) {
        q <- c(1L, 1L, 2L, 2L)[w > 0]
        s <- c(1L, 2L, 1L, 2L)[w > 0]
        w <- w[w > 0]
        expect_identical(SegStability:::.geometryMatching(q, s, w),
                         wp09_oracle(q, s, w))
    }
    expect_identical(SegStability:::.geometryMatching(integer(), integer(),
                                                    numeric()), integer())
    # Both scalar totals round to 1, but the diagonal has a real improvement.
    expect_identical(sum(c(1, 2^-54)), sum(c(.5, .5)))
    expect_identical(SegStability:::.geometryMatching(c(1L, 1L, 2L, 2L),
        c(1L, 2L, 1L, 2L), c(1, .5, .5, 2^-54)), c(1L, 4L))
    # Reversing that improvement must override the lexicographically first pair.
    expect_identical(SegStability:::.geometryMatching(c(1L, 1L, 2L, 2L),
        c(1L, 2L, 1L, 2L), c(.5, 1, 2^-54, .5)), c(2L, 3L))
})

test_that("tiny exhaustive oracle covers sparse rectangular matchings", {
    for (mask in 0:63) {
        selected <- as.logical(intToBits(mask)[1:6])
        q <- rep(1:2, each = 3)[selected]
        s <- rep(1:3, 2)[selected]
        for (weights in list(rep(.25, 6), c(.5, .25, .125, .75, .5, .25))) {
            w <- weights[selected]
            expect_identical(SegStability:::.geometryMatching(q, s, w),
                             wp09_oracle(q, s, w))
        }
    }
})
