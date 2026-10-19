# Assisted-by: OpenAI Codex
.provenanceFields <- function() {
    c("source", "source_version", "retrieval_date", "source_hashes",
      "software", "software_version", "software_hashes", "command",
      "parameters", "model", "threads", "runtime", "observation_window",
      "output_hashes")
}

.structuredValue <- function(x) {
    is.list(x) || methods::is(x, "SimpleList") || methods::is(x, "DataFrame")
}

.unknownMetadata <- function(x) {
    if (is.null(x) || !length(x)) {
        return(TRUE)
    }
    if (.structuredValue(x)) {
        return(any(vapply(as.list(x), .unknownMetadata, logical(1))))
    }
    if (is.atomic(x)) {
        return(anyNA(x))
    }
    FALSE
}

.validationWindow <- function(x) {
    .baseMeasurement(x) && length(x) == 4L &&
        identical(sort(names(x)), c("xmax", "xmin", "ymax", "ymin")) &&
        all(is.finite(x)) && x[["xmin"]] < x[["xmax"]] &&
        x[["ymin"]] < x[["ymax"]]
}

.hashMetadata <- function(x) {
    if (is.character(x) && !is.object(x)) {
        return(all(is.na(x) | nzchar(x)))
    }
    if (!(is.list(x) || methods::is(x, "SimpleList"))) {
        return(FALSE)
    }
    if (!length(x)) {
        return(TRUE)
    }
    nm <- names(x)
    !is.null(nm) && !anyNA(nm) && all(nzchar(nm)) &&
        !anyDuplicated(nm) && all(vapply(as.list(x), .hashMetadata, logical(1)))
}

.provenanceShape <- function(field, x) {
    if (field %in% c("source_hashes", "software_hashes", "output_hashes")) {
        return(.hashMetadata(x))
    }
    if (field %in% c("parameters", "model")) {
        return((is.character(x) && !is.object(x)) || .structuredValue(x))
    }
    if (field == "threads") {
        return(identical(class(x), "integer") && length(x) == 1L &&
                   (is.na(x) || x > 0L))
    }
    if (field == "runtime") {
        return(.baseMeasurement(x) && length(x) == 1L && !is.nan(x) &&
                   (is.na(x) || (is.finite(x) && x >= 0)))
    }
    if (field == "observation_window") {
        return((methods::is(x, "SimpleList") && !length(x)) ||
                   .validationWindow(x))
    }
    is.character(x) && !is.object(x) && length(x) == 1L &&
        (is.na(x) || nzchar(x))
}

.provenanceChecks <- function(object) {
    if (methods::is(object, "SegmentationPartition") ||
        methods::is(object, "SegmentationSet") ||
        methods::is(object, "SegmentationComparison")) {
        object <- provenance(object)
    }
    structure <- .validationCheck("provenance.structure",
        if (.namedValidationList(object)) TRUE else
            "provenance must be a uniquely named SimpleList")
    out <- list(structure)
    if (.validationOK(structure)) {
        for (field in .provenanceFields()) {
            check <- paste0("provenance.", field)
            if (!field %in% names(object)) {
                row <- .validationReport(check, "FAIL", "error",
                                         "Required provenance key is missing")
            } else if (!.provenanceShape(field, object[[field]])) {
                row <- .validationReport(check, "FAIL", "error",
                                         "Malformed provenance field value")
            } else if (.unknownMetadata(object[[field]])) {
                row <- .validationUnknown(check, "Required evidence is unknown")
            } else {
                row <- .validationReport(check, "PASS",
                    message = "Field is present and has the expected shape")
            }
            out[[length(out) + 1L]] <- row
        }
    } else {
        out[[length(out) + 1L]] <- .validationUnknown("provenance.fields",
            "Field checks unavailable because provenance structure failed")
    }
    out[[length(out) + 1L]] <- .validationUnknown("provenance.authenticity",
        "Source authenticity, source correctness and truth are not attested")
    .validationBind(out)
}

#' Validate provenance shape and completeness
#'
#' Checks standardized fields without completing or changing supplied metadata.
#' Unknown evidence remains unknown; completeness does not attest authenticity.
#' @param object A SimpleList or a SegmentationPartition, SegmentationSet or
#'   SegmentationComparison, whose provenance accessor supplies the metadata.
#' @param strict Logical scalar. Stop only if a report contains FAIL/error.
#'   UNKNOWN and FAIL/warning findings do not themselves cause an error.
#' @return A DataFrame with character columns check, status, severity and
#'   message.
#'   Status is PASS, FAIL, UNKNOWN or NOT_APPLICABLE; severity is error, warning
#'   or note. Missing required keys and malformed values are FAIL/error;
#'   typed unknowns and appropriate empty structured values are UNKNOWN/note.
#' @details
#' Required keys are source, source_version, retrieval_date, source_hashes,
#' software, software_version, software_hashes, command, parameters, model,
#' threads, runtime, observation_window and output_hashes. Text fields are
#' character scalars; hashes are character vectors or uniquely named nested
#' list/SimpleList records. Parameters/model accept character or structured
#' values. Threads is a positive base integer scalar; runtime is a finite
#' nonnegative base numeric scalar. Typed NA values remain unknown. A window
#' is a finite named numeric xmin/xmax/ymin/ymax bbox with increasing bounds,
#' or an empty SimpleList unknown placeholder. Extra keys are retained.
#' Reports are deterministic and side-effect free. Rows describe checks rather
#' than an aggregate validity or authenticity flag. FAIL/error closes the
#' relevant gate; UNKNOWN indicates unavailable evidence. Strict mode throws
#' only for FAIL/error, never solely for FAIL/warning or UNKNOWN. Backend
#' warnings are captured in diagnostics; a warned check without a failure is
#' UNKNOWN/warning rather than an unsupported PASS.
#' @examples
#' p <- SegmentationPartition(
#'     S4Vectors::DataFrame(molecule_id = character(),
#'         feature_id = character()),
#'     "empty", S4Vectors::DataFrame(cell_id = character()),
#'     character(), character())
#' validateProvenance(p, strict = TRUE)
#' @export
validateProvenance <- function(object, strict = FALSE) {
    .validationFinish(.provenanceChecks(object), strict)
}
