.coordinateShape <- function(field, x) {
    if (field %in% c("units", "specimen_id", "crs")) {
        return(is.character(x) && !is.object(x) && length(x) == 1L &&
                   (is.na(x) || nzchar(x)))
    }
    if (field %in% c("origin", "axis_directions")) {
        return((is.character(x) && !is.object(x)) ||
                   (.baseMeasurement(x) && all(is.finite(x) | is.na(x)) &&
                    !any(is.nan(x))))
    }
    (is.character(x) && !is.object(x)) || .structuredValue(x)
}

.coordinateListProblems <- function(info) {
    if (!.namedValidationList(info)) {
        return("coordinate_info must be a uniquely named SimpleList")
    }
    fields <- c("units", "origin", "axis_directions", "transforms",
                "registration", "specimen_id", "crs")
    for (field in intersect(fields, names(info))) {
        if (!.coordinateShape(field, info[[field]])) {
            return(paste("Malformed coordinate declaration:", field))
        }
    }
    TRUE
}

.coordinateDeclarations <- function(info, stored = NULL, reference = NULL,
                                    specimen = NULL, prefix = "coordinates") {
    if (is.null(info)) {
        info <- if (is.null(stored)) S4Vectors::SimpleList() else stored
    }
    structure <- .validationCheck(paste0(prefix, ".declarations"),
        if (.namedValidationList(info)) TRUE else
            "coordinate_info must be a uniquely named SimpleList")
    out <- list(structure)
    for (name in c("stored", "reference")) {
        candidate <- if (name == "stored") stored else reference
        if (!is.null(candidate)) {
            out[[length(out) + 1L]] <- .validationCheck(
                paste0(prefix, ".declarations.", name),
                .coordinateListProblems(candidate))
        }
    }
    fields <- c("units", "origin", "axis_directions", "transforms",
                "registration", "specimen_id", "crs")
    if (.validationOK(structure)) {
        for (field in fields) {
            value <- info[[field]]
            if (field == "specimen_id" && is.null(value) &&
                !is.null(specimen)) {
                value <- specimen
            }
            check <- paste0(prefix, ".", field)
            known_none <- field == "transforms" &&
                .structuredValue(value) && !length(value)
            row <- if (is.null(value)) {
                .validationUnknown(check, "Coordinate declaration is absent")
            } else if (!.coordinateShape(field, value)) {
                .validationReport(check, "FAIL", "error",
                                  "Malformed coordinate declaration")
            } else if (.unknownMetadata(value) && !known_none) {
                .validationUnknown(check, "Coordinate declaration is unknown")
            } else {
                .validationReport(check, "PASS",
                    message = "Declaration shape checked; source not attested")
            }
            out[[length(out) + 1L]] <- row
            candidates <- list(stored = stored, reference = reference)
            for (name in names(candidates)) {
                other <- candidates[[name]]
                if (.namedValidationList(other) && field %in% names(other) &&
                    !is.null(value) && .coordinateShape(field, value) &&
                    (!.unknownMetadata(value) || known_none)) {
                    expected <- other[[field]]
                    expected_none <- field == "transforms" &&
                        .structuredValue(expected) && !length(expected)
                    if (!.unknownMetadata(expected) || expected_none) {
                        out[[length(out) + 1L]] <- .validationCheck(
                            paste0(check, ".", name),
                            if (identical(value, expected)) TRUE else
                                "Known coordinate declarations conflict")
                    }
                }
            }
            if (field == "specimen_id" && !is.null(specimen) &&
                !is.null(value) && !.unknownMetadata(value)) {
                out[[length(out) + 1L]] <- .validationCheck(
                    paste0(check, ".owner"),
                    if (identical(value, specimen)) TRUE else
                        "Specimen declaration conflicts with owning set")
            }
        }
    }
    out[[length(out) + 1L]] <- .validationUnknown(
        paste0(prefix, ".registration_evidence"),
        "Declarations alone do not verify registration or source linkage")
    .validationBind(out)
}

.coordinateReference <- function(reference) {
    .validationCheck("coordinates.reference", {
        if (!methods::is(reference, "SegmentationSet")) {
            stop("reference must be a SegmentationSet")
        }
        valid <- methods::validObject(reference, test = TRUE, complete = TRUE)
        if (!isTRUE(valid)) {
            stop(paste(valid, collapse = "; "))
        }
        declarations <- .coordinateListProblems(reference@coordinate_info)
        if (!isTRUE(declarations)) {
            stop(declarations)
        }
        declared <- reference@coordinate_info$specimen_id
        if (!is.null(declared) && !.unknownMetadata(declared) &&
            !identical(declared, reference@specimen_id)) {
            stop("Reference specimen declaration conflicts with owning set")
        }
        TRUE
    })
}

#' Validate molecule coordinates without transformation
#'
#' Finite base numeric/double or integer spatial measurements are accepted.
#' The no-double identifier rule applies to molecule/cell IDs, not coordinates.
#' Exact ID alignment may index temporary diagnostics, but never reorders input.
#' @param object A SegmentationSet or DataFrame with molecule_id, x, y and
#'   optional z columns. NULL denotes absent optional coordinates.
#' @param coordinate_info Optional uniquely named SimpleList declaring units,
#'   origin, axis_directions, transforms, registration, specimen_id and crs.
#'   Missing declarations remain unknown. An empty structured transforms value
#'   explicitly means none. Known conflicting declarations fail; units and
#'   declarations are compared exactly, without aliases or inferred transforms.
#' @param reference Optional SegmentationSet for exact molecule coverage,
#'   coordinate values and declaration consistency checks.
#' @param strict Logical scalar. Stop only for FAIL/error findings.
#' @return A DataFrame with character check, status, severity and message
#'   columns.
#'   Absent coordinates are NOT_APPLICABLE/note. Malformed supplied coordinates
#'   are FAIL/error. Matching values/declarations do not attest registration.
#' @details
#' Units, specimen_id and crs are character scalars. Origin and axis_directions
#' accept character or base numeric vectors. Transforms/registration accept
#' character or structured records; their contents are not executed or attested.
#' Explicit declarations must agree with known stored/reference declarations.
#' Finite measurements are checked without inferred range limits or imputation.
#' Source identity linkage, registration verification and authenticity remain
#' unchecked even when values and declarations agree.
#' @examples
#' d <- S4Vectors::DataFrame(molecule_id = c("m1", "m2"),
#'     x = c(0.25, 1.5), y = c(2.5, 3.75))
#' validateCoordinates(d, strict = TRUE)
#' @export
validateCoordinates <- function(object, coordinate_info = NULL,
                                reference = NULL, strict = FALSE) {
    if (is.null(object)) {
        return(.validationFinish(.validationReport("coordinates.presence",
            "NOT_APPLICABLE", message = "Optional coordinates not supplied"),
            strict))
    }
    is_set <- methods::is(object, "SegmentationSet")
    data <- if (is_set) object@molecules else object
    input <- .validationCheck("coordinates.input",
        if (methods::is(data, "DataFrame")) TRUE else
            "Expected a SegmentationSet or DataFrame")
    if (!.validationOK(input)) {
        return(.validationFinish(input, strict))
    }
    out <- list(input)
    present <- intersect(c("x", "y", "z"), colnames(data))
    if (!length(present)) {
        out[[length(out) + 1L]] <- .validationReport("coordinates.presence",
            "NOT_APPLICABLE", message = "Optional coordinates are absent")
    } else {
        shape <- .validationCheck("coordinates.shape", {
            if (!all(c("molecule_id", "x", "y") %in% colnames(data)) ||
                anyDuplicated(colnames(data))) {
                stop("Coordinates require unique molecule_id, x and y columns")
            }
            TRUE
        })
        identity <- .validationCheck("coordinates.identity", {
            .validateIds(data$molecule_id, "coordinate molecule_id")
            TRUE
        })
        out <- c(out, list(shape, identity))
        for (axis in present) {
            out[[length(out) + 1L]] <- .validationCheck(
                paste0("coordinates.values.", axis), {
                    v <- data[[axis]]
                    if (!.baseMeasurement(v) || !is.null(dim(v)) ||
                        length(v) != nrow(data) || any(!is.finite(v))) {
                        stop("Coordinates require finite base numeric vectors")
                    }
                    TRUE
                })
        }
        if (!is.null(reference)) {
            ref <- .coordinateReference(reference)
            out[[length(out) + 1L]] <- ref
            if (.validationOK(ref) && .validationOK(identity) &&
                .validationOK(shape)) {
                out[[length(out) + 1L]] <- .validationCheck(
                    "coordinates.alignment", {
                        idx <- .alignMolecules(data$molecule_id,
                                               reference@molecules$molecule_id)
                        for (axis in present) {
                            other <- reference@molecules[[axis]]
                            if (is.null(other)) {
                                stop("Reference coordinate dimension missing")
                            }
                            if (!.baseMeasurement(other) ||
                                any(!is.finite(other)) ||
                                !all(data[[axis]][idx] == other)) {
                                stop("Identity-aligned coordinates differ")
                            }
                        }
                        if (!identical(present, intersect(c("x", "y", "z"),
                            colnames(reference@molecules)))) {
                            stop("Reference coordinate dimensions differ")
                        }
                        TRUE
                    })
            } else {
                out[[length(out) + 1L]] <- .validationUnknown(
                    "coordinates.alignment",
                    "Reference/alignment checks failed")
            }
        }
    }
    stored <- if (is_set) object@coordinate_info else NULL
    specimen <- if (is_set) object@specimen_id else NULL
    ref_info <- NULL
    if (methods::is(reference, "SegmentationSet")) {
        ref_info <- reference@coordinate_info
        ref_info$specimen_id <- reference@specimen_id
    }
    out[[length(out) + 1L]] <- .coordinateDeclarations(coordinate_info,
        stored, ref_info, specimen)
    .validationFinish(.validationBind(out), strict)
}
