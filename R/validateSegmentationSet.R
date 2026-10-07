#' Validate a shared-universe segmentation set explicitly
#'
#' Checks exact identity, feature mapping, full partition binding and
#'   conservation.
#' Recomputes the shared universe once and each assignment fingerprint once when
#' enabled. Coordinates are optional; this function checks declaration
#'   structure,
#' without running spatial operations. No comparisons or repairs are performed.
#' @param object A SegmentationSet.
#' @param check_integrity Logical scalar. Recompute existing fingerprints;
#'   FALSE reports integrity as UNKNOWN rather than claiming it was checked.
#' @param strict Logical scalar. Stop only for FAIL/error findings.
#' @return A DataFrame with character check, status, severity and message
#'   columns.
#'   Nested findings follow partition order. UNKNOWN evidence cannot certify a
#'   gate. Integrity fingerprints do not attest source authenticity or
#'   provenance.
#' @examples
#' m <- S4Vectors::DataFrame(molecule_id = "m1", feature_id = "g1")
#' p <- SegmentationPartition(m, "s", S4Vectors::DataFrame(cell_id = "c1"),
#'     "m1", "assigned", "c1")
#' x <- SegmentationSet(m, list(s = p), "specimen")
#' validateSegmentationSet(x, strict = TRUE)
#' @export
validateSegmentationSet <- function(object, check_integrity = TRUE,
                                    strict = FALSE) {
    .validationFlag(check_integrity)
    input <- .validationCheck("set.input",
        if (methods::is(object, "SegmentationSet")) TRUE else
            "Expected a SegmentationSet")
    if (!.validationOK(input)) {
        return(.validationFinish(input, strict))
    }
    structure <- .validationCheck("set.structure",
        methods::validObject(object, test = TRUE, complete = TRUE))
    universe <- .validationCheck("universe.structure", {
        .checkUniverseStructure(object@molecules, object@features)
        TRUE
    })
    specimen <- .validationCheck("set.specimen_id",
        if (.scalarString(object@specimen_id)) TRUE else
            "specimen_id must be one nonmissing nonempty character string")
    out <- list(input, structure, universe, specimen)
    out[[length(out) + 1L]] <- if (check_integrity &&
        .validationOK(universe)) {
        .validationCheck("integrity.universe", {
            .validateUniverse(S4Vectors::SimpleList(
                molecules = object@molecules, features = object@features,
                fingerprint_value = object@universe_fingerprint,
                fingerprint_version = object@fingerprint_version,
                fingerprint_algorithm = object@fingerprint_algorithm))
            TRUE
        })
    } else {
        .validationUnknown("integrity.universe",
            "Recomputation disabled or universe structure failed")
    }
    for (i in seq_along(object@segmentations)) {
        p <- object@segmentations[[i]]
        prefix <- paste0("partition[", i, "].")
        binding <- .validationCheck("universe.binding",
            .partitionBindingProblems(p, nrow(object@molecules),
                object@universe_fingerprint, object@fingerprint_version,
                object@fingerprint_algorithm))
        report <- .partitionValidation(p, check_integrity)
        if (.validationOK(binding)) {
            row <- report$check == "universe.verification"
            report$status[row] <- "PASS"
            report$message[row] <- paste("Stored binding and length agree with",
                                         "owning set; see integrity.universe")
        }
        out[[length(out) + 1L]] <- .validationPrefix(
            .validationBind(list(binding, report)), prefix)
    }
    for (field in c("coordinate_info", "metadata")) {
        value <- if (field == "coordinate_info") {
            .coordinateListProblems(object@coordinate_info)
        } else if (methods::is(object@metadata, "SimpleList")) {
            TRUE
        } else {
            "metadata must be a SimpleList"
        }
        out[[length(out) + 1L]] <- .validationCheck(
            paste0("set.", field), value)
    }
    out[[length(out) + 1L]] <- .provenanceChecks(object)
    .validationFinish(.validationBind(out), strict)
}
