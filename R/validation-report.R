# Assisted-by: OpenAI Codex
.validationReport <- function(check, status, severity = "note", message) {
    S4Vectors::DataFrame(check = check, status = status, severity = severity,
                        message = message)
}

.validationCheck <- function(check, expr, severity = "error") {
    evidence <- new.env(parent = emptyenv())
    evidence$warnings <- character()
    result <- tryCatch(withCallingHandlers({
        value <- expr
        if (!isTRUE(value)) {
            stop(paste(value, collapse = "; "), call. = FALSE)
        }
        NULL
    }, warning = function(w) {
        evidence$warnings <- c(evidence$warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
    if (!is.null(result)) {
        .validationReport(check, "FAIL", severity,
            paste(c(result, evidence$warnings), collapse = "; "))
    } else if (length(evidence$warnings)) {
        .validationReport(check, "UNKNOWN", "warning",
            paste(evidence$warnings, collapse = "; "))
    } else {
        .validationReport(check, "PASS", message = "Check passed")
    }
}

.validationBind <- function(reports) {
    do.call(rbind, reports)
}

.validationOK <- function(report) all(report$status == "PASS")

.validationFinish <- function(report, strict) {
    if (!identical(strict, TRUE) && !identical(strict, FALSE)) {
        stop("strict must be TRUE or FALSE", call. = FALSE)
    }
    failed <- report$status == "FAIL" & report$severity == "error"
    if (strict && any(failed)) {
        stop("Validation failed: ", paste(report$check[failed],
                                         collapse = ", "), call. = FALSE)
    }
    report
}

.validationFlag <- function(x) {
    if (!identical(x, TRUE) && !identical(x, FALSE)) {
        stop("check_integrity must be TRUE or FALSE", call. = FALSE)
    }
}

.namedValidationList <- function(x) {
    methods::is(x, "SimpleList") && (!length(x) ||
        (!is.null(names(x)) && !anyNA(names(x)) &&
         all(nzchar(names(x))) && !anyDuplicated(names(x))))
}

.validationUnknown <- function(check, message) {
    .validationReport(check, "UNKNOWN", message = message)
}

.validationPrefix <- function(report, prefix) {
    report$check <- paste0(prefix, report$check)
    report
}

.baseMeasurement <- function(x) {
    (identical(class(x), "numeric") && typeof(x) == "double") ||
        (identical(class(x), "integer") && typeof(x) == "integer")
}
