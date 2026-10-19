## Scientific IDs retain their family; keys are only for certified exact paths.
.idFamily <- function(x, label = "ID") {
    kind <- as.character(class(x))
    if (identical(kind, "integer") && typeof(x) == "integer") {
        return("integer")
    }
    if (identical(kind, "character") && typeof(x) == "character") {
        return("character")
    }
    if (identical(kind, "LLint") && methods::is(x, "LLint")) {
        return("LLint")
    }
    if (identical(kind, "integer64") && typeof(x) == "double") {
        if (!requireNamespace("bit64", quietly = TRUE)) {
            stop(label, ": install the optional 'bit64' package to use integer64 IDs",
                 call. = FALSE)
        }
        return("integer64")
    }
    stop(label, ": unsupported ID representation", call. = FALSE)
}

.idKeys <- function(x) {
    .idFamily(x)
    ## Both 64-bit decimal conversions were independently audited at Gate D.
    unname(as.character(x))
}

.validateIds <- function(x, label = "ID", unique = TRUE) {
    kind <- .idFamily(x, label)
    if (any(is.na(x))) {
        stop(label, ": Missing IDs", call. = FALSE)
    }
    if (kind == "character" && any(!nzchar(x))) {
        stop(label, ": Empty IDs", call. = FALSE)
    }
    if (unique && anyDuplicated(.idKeys(x))) {
        stop(label, ": Duplicate IDs", call. = FALSE)
    }
    x
}

.matchIds <- function(x, table) {
    .validateIds(x, unique = FALSE)
    .validateIds(table)
    if (.idFamily(x) != .idFamily(table)) {
        stop("ID families must agree; mixed-family matching is unsupported",
             call. = FALSE)
    }
    base::match(.idKeys(x), .idKeys(table))
}

.idOrder <- function(x) {
    .validateIds(x, unique = FALSE)
    switch(.idFamily(x),
        integer = base::order(x, method = "radix"),
        character = base::order(x, method = "radix"),
        LLint = base::order(x),
        integer64 = bit64::order(x))
}

.alignMolecules <- function(external, universe) {
    .validateIds(external, "assignment molecule_id")
    permutation <- .matchIds(universe, external)
    extra <- is.na(.matchIds(external, universe))
    if (anyNA(permutation) || any(extra)) {
        stop("Assignment universe mismatch: ", sum(is.na(permutation)),
             " missing, ", sum(extra), " extra molecules", call. = FALSE)
    }
    permutation
}

.validateFeatures <- function(x, unique = FALSE) {
    if (!identical(class(x), "character") || typeof(x) != "character") {
        stop("feature_id must be base character; no automatic coercion",
             call. = FALSE)
    }
    if (anyNA(x) || any(!nzchar(x))) {
        stop("feature_id must not contain missing or empty identities",
             call. = FALSE)
    }
    if (unique && anyDuplicated(x)) {
        stop("Duplicate feature_id dictionary entries", call. = FALSE)
    }
    x
}

.checkIndexSize <- function(n) {
    if (n > .Machine$integer.max) {
        stop("Population exceeds the supported integer-index range",
             call. = FALSE)
    }
    invisible(NULL)
}
