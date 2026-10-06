for (kind in c("integer", "character", "LLint", "integer64")) {
    local({
        representation <- kind
        test_that(paste(representation, "preserves every laboratory identity path"), {
            if (representation == "integer64") {
                skip_if_not_installed("bit64")
            }
            fixture <- id_lab_fixture(representation)
            x <- fixture$x
            expect_id_lab_exact(id_lab_validate(x), x)
            expect_identical(as.character(x), fixture$values)
            expect_true(all(x == x))
            expect_false(any(x[1L] == x[2L]))

            duplicate <- x[c(1L, 2L, 1L)]
            keys <- if (representation == "LLint") {
                as.character(duplicate)
            } else {
                duplicate
            }
            expect_identical(base::duplicated(keys), c(FALSE, FALSE, TRUE))
            expect_error(id_lab_validate(duplicate), "Duplicate IDs")
            missing <- if (representation == "LLint") {
                # LLint does not support missing subscripts; test missing values directly.
                S4Vectors::as.LLint(c(fixture$values[1L], NA_character_))
            } else {
                x[c(1L, NA_integer_)]
            }
            expect_error(id_lab_validate(missing), "Missing IDs")

            permutation <- rev(seq_along(x))
            query <- x[permutation]
            expect_identical(base::match(query, x), permutation)
            expect_identical(base::match(x[1:2], x[-1L]), c(NA_integer_, 1L))
            expect_identical(base::match(query, x)[permutation], seq_along(x))
            expect_id_lab_exact(query[permutation], x)

            ordering <- if (representation == "integer64") {
                bit64::order(x)
            } else if (representation == "character") {
                base::order(x, method = "radix")
            } else {
                base::order(x)
            }
            expect_identical(ordering, fixture$sorted)
            expect_id_lab_exact(x[ordering], x[fixture$sorted])
            sorted <- if (representation == "integer64") {
                base::sort(x)
            } else if (representation == "character") {
                base::sort(x, method = "radix")
            } else {
                base::sort(x)
            }
            expect_id_lab_exact(sorted, x[fixture$sorted])

            frame <- S4Vectors::DataFrame(id = x, position = seq_along(x))
            expect_id_lab_exact(frame$id, x)
            expect_id_lab_exact(frame[["id"]], x)
            reordered <- frame[permutation, , drop = FALSE]
            expect_id_lab_exact(reordered$id, query)
            aligned <- reordered[base::match(x, reordered$id), , drop = FALSE]
            expect_id_lab_exact(aligned$id, x)
            expect_identical(aligned$position, seq_along(x))

            path <- tempfile(fileext = ".rds")
            on.exit(unlink(path), add = TRUE)
            for (object in list(x, frame, reordered)) {
                saveRDS(object, path)
                restored <- readRDS(path)
                expect_true(identical(restored, object,
                                      num.eq = FALSE, single.NA = FALSE))
                ids <- if (methods::is(restored, "DataFrame")) restored$id else restored
                expected <- if (methods::is(object, "DataFrame")) object$id else object
                expect_id_lab_exact(ids, expected)
            }
        })
    })
}

test_that("LLint native duplication is unsupported but decimal keys round-trip", {
    x <- id_lab_fixture("LLint")$x
    expect_error(base::duplicated(x), "applies only to vectors")
    expect_error(x[NA_integer_], "NA")
    expect_id_lab_exact(S4Vectors::as.LLint(as.character(x)), x)
    for (value in c("-9223372036854775808", "9223372036854775808")) {
        expect_warning(outside <- S4Vectors::as.LLint(value), "out-of-range")
        expect_error(id_lab_validate(outside), "Missing IDs")
    }
})

test_that("integer64 needs explicit ordering and exact assertions", {
    skip_if_not_installed("bit64")
    fixture <- id_lab_fixture("integer64")
    x <- fixture$x
    expect_false(identical(base::order(x), fixture$sorted))
    expect_identical(bit64::order(x), fixture$sorted)
    expect_id_lab_exact(bit64::as.integer64(as.character(x)), x)
    # Default relative tolerance accepts these distinct scientific identities.
    expect_success(expect_equal(x[1L], x[2L]))
    expect_failure(expect_equal(x[1L], x[2L], tolerance = 0))
    expect_failure(expect_identical(x[1L], x[2L]))
    expect_failure(expect_identical(as.character(x[1L]), as.character(x[2L])))
    for (value in c("-9223372036854775808", "9223372036854775808")) {
        expect_warning(outside <- bit64::as.integer64(value), "integer64 range")
        expect_error(id_lab_validate(outside), "Missing IDs")
    }
})

test_that("LLint exact assertions distinguish adjacent values", {
    x <- id_lab_fixture("LLint")$x
    expect_failure(expect_equal(x[1L], x[2L]))
    expect_failure(expect_identical(x[1L], x[2L]))
})

test_that("audited 64-bit paths never call numeric conversion methods", {
    ll_namespace <- asNamespace("S4Vectors")
    trace(".from_LLint_to_numeric", where = ll_namespace,
          tracer = quote(stop("Forbidden numeric ID conversion")), print = FALSE)
    on.exit(untrace(".from_LLint_to_numeric", where = ll_namespace), add = TRUE)
    representations <- "LLint"
    if (requireNamespace("bit64", quietly = TRUE)) {
        bit_namespace <- asNamespace("bit64")
        trace("as.double.integer64", where = bit_namespace,
              tracer = quote(stop("Forbidden numeric ID conversion")), print = FALSE)
        on.exit(untrace("as.double.integer64", where = bit_namespace), add = TRUE)
        representations <- c(representations, "integer64")
    }
    for (representation in representations) {
        x <- id_lab_fixture(representation)$x
        expect_id_lab_exact(id_lab_validate(x), x)
        expect_identical(base::match(x[2:1], x), c(2L, 1L))
        expect_true(all(x == x))
        ordering <- if (representation == "LLint") base::order(x) else bit64::order(x)
        expect_id_lab_exact(x[ordering], base::sort(x))
        frame <- S4Vectors::DataFrame(id = x)
        expect_id_lab_exact(frame[c(2L, 1L), , drop = FALSE]$id, x[2:1])
        path <- tempfile(fileext = ".rds")
        on.exit(unlink(path), add = TRUE)
        saveRDS(x, path)
        expect_id_lab_exact(readRDS(path), x)
    }
})
