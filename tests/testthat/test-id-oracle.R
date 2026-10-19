# Assisted-by: OpenAI Codex
# Identity cases ported from scripts/144_validate_final_portability_cases.py
# (key(), lines 20-23; __main__, lines 73-76 and 84-85). No graph logic is ported.
test_that("historical integer, large, character and UUID identities remain exact", {
    fixtures <- list(
        c(1L, 2L),
        c("9223372036854775809", "9223372036854775810"),
        c("a", "b"),
        c("550e8400-e29b-41d4-a716-446655440000",
          "550e8400-e29b-41d4-a716-446655440001")
    )
    for (x in fixtures) {
        expect_id_lab_exact(id_lab_validate(x), x)
        expect_identical(base::match(x[2:1], x), c(2L, 1L))
        expect_id_lab_exact(x[2:1][2:1], x)
    }
    # The historical mixed-type case represents separate typed identities in R.
    expect_false(identical(1L, "1"))
    expect_false(identical(list(kind = "integer", id = 1L),
                           list(kind = "character", id = "1")))
    # Base matching across families conflates spellings; it is not a supported path.
    expect_identical(base::match(1L, "1"), 1L)
    expect_error(id_lab_validate(list(1L, "1")), "Unsupported")
})

test_that("historical unsafe, missing and duplicate identities fail", {
    # Includes the oracle's float(2**63), None, True and empty string cases.
    rejected <- list(2^63, NA_integer_, TRUE, "", 1, c(1, 2),
        NA_character_, NA_real_, NaN, Inf, -Inf, factor("a"),
        1 + 1i, as.raw(1L), structure(1L, class = "unknown_id"), NULL)
    for (x in rejected) {
        expect_error(id_lab_validate(x))
    }
    # Duplicate census and transcript IDs have the same uniqueness invariant.
    expect_error(id_lab_validate(c(1L, 1L)), "Duplicate IDs")
    expect_error(id_lab_validate(c(7L, 7L)), "Duplicate IDs")
    expect_error(id_lab_validate(c("a", "a")), "Duplicate IDs")
    expect_id_lab_exact(id_lab_validate(c("001", "1", " 1 ", "1e0")),
                        c("001", "1", " 1 ", "1e0"))
})
