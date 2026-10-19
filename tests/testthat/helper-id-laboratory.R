# Assisted-by: OpenAI Codex
# Test-only acceptance probe; this is not a production ID utility or API.
id_lab_validate <- function(x) {
    kind <- as.character(class(x))
    allowed <- identical(kind, "integer") || identical(kind, "character") ||
        (identical(kind, "LLint") && isS4(x)) ||
        (identical(kind, "integer64") && typeof(x) == "double")
    if (!allowed) {
        stop("Unsupported ID representation")
    }
    if (any(is.na(x))) {
        stop("Missing IDs")
    }
    if (identical(kind, "character") && any(!nzchar(x))) {
        stop("Empty IDs")
    }
    # LLint lacks duplicated(); its exact decimal conversion is audited separately.
    keys <- if (identical(kind, "LLint")) as.character(x) else x
    if (any(base::duplicated(keys))) {
        stop("Duplicate IDs")
    }
    invisible(x)
}

id_lab_large_values <- function() {
    c("9007199254740993", "9007199254740992", "9223372036854775807",
      "-9223372036854775807", "0", "9007199254740994",
      "9223372036854775806", "-1", "-9223372036854775806")
}

id_lab_fixture <- function(kind) {
    values <- switch(kind,
        integer = c("2", "2147483647", "-2147483647", "0"),
        character = c("b", "a",
            "550e8400-e29b-41d4-a716-446655440000",
            "550e8400-e29b-41d4-a716-446655440001",
            "9223372036854775809", "9223372036854775810", "001", "1"),
        id_lab_large_values()
    )
    x <- switch(kind,
        integer = as.integer(values),
        character = values,
        LLint = S4Vectors::as.LLint(values),
        integer64 = bit64::as.integer64(values)
    )
    sorted <- switch(kind,
        integer = c(3L, 4L, 1L, 2L),
        character = c(7L, 8L, 3L, 4L, 5L, 6L, 2L, 1L),
        c(4L, 9L, 8L, 5L, 2L, 1L, 6L, 7L, 3L)
    )
    list(x = x, values = values, sorted = sorted)
}

# Character values and bit-exact identity avoid numeric assertion tolerance.
# Serialization format itself can vary with ALTREP storage of identical integers.
expect_id_lab_exact <- function(actual, expected) {
    expect_identical(class(actual), class(expected))
    expect_identical(as.character(actual), as.character(expected))
    expect_true(identical(actual, expected, num.eq = FALSE, single.NA = FALSE))
}
