.partitionValidation <- function(object, check_integrity) {
    input <- .validationCheck("segmentation.input",
        if (methods::is(object, "SegmentationPartition")) TRUE else
            "Expected a SegmentationPartition")
    if (!.validationOK(input)) {
        return(input)
    }
    structure <- .validationCheck("segmentation.structure",
        methods::validObject(object, test = TRUE, complete = TRUE))
    census <- .validationCheck("census.identity", {
        if (!"cell_id" %in% colnames(object@cells) ||
            anyDuplicated(colnames(object@cells))) {
            stop("Census requires an unambiguous cell_id column")
        }
        .validateIds(object@cells$cell_id, "cell census")
        TRUE
    })
    status <- .validationCheck("assignment.status", {
        s <- object@status
        if (length(s) != length(object@cell_index) || anyNA(s) ||
            any(!s %in% c(1L, 2L, 3L))) {
            stop("Invalid status codes or assignment lengths")
        }
        TRUE
    })
    out <- list(input, structure, census, status,
        .validationUnknown("census.source_completeness",
            "External census completeness cannot be inferred from assignments"))
    if (.validationOK(status) && .validationOK(census)) {
        assigned <- .validationCheck("assignment.cells", {
            a <- object@status == 1L
            i <- object@cell_index
            if (anyNA(i[a]) || any(i[a] < 1L | i[a] > nrow(object@cells)) ||
                any(!is.na(i[!a]))) {
                stop("Invalid cell references or forced unresolved assignment")
            }
            TRUE
        })
        conservation <- .validationCheck("assignment.conservation", {
            counts <- tabulate(object@status, nbins = 3L)
            mass <- tabulate(object@cell_index[object@status == 1L],
                             nbins = nrow(object@cells))
            if (sum(counts) != length(object@status) ||
                sum(mass) != counts[1L]) {
                stop("Assignment/status conservation failed")
            }
            TRUE
        })
        ambiguity <- .validationCheck("ambiguity.consistency",
            .ambiguityProblems(object@ambiguity, object@status,
                               nrow(object@cells)))
        out <- c(out, list(assigned, conservation, ambiguity))
        if (.validationOK(ambiguity)) {
            evidence <- object@ambiguity$evidence
            incomplete <- any(!evidence$candidate_ids_known) ||
                anyNA(evidence$candidate_count)
            out[[length(out) + 1L]] <- if (incomplete) {
                .validationUnknown("ambiguity.completeness",
                    "Candidate identities or counts remain unknown")
            } else {
                .validationReport("ambiguity.completeness", "PASS",
                                  message = "Candidate evidence is complete")
            }
        }
    } else {
        out[[length(out) + 1L]] <- .validationUnknown("assignment.scientific",
            "Assignment checks unavailable after status/census failure")
    }
    out[[length(out) + 1L]] <- if (check_integrity &&
        .validationOK(structure) && .validationOK(census)) {
        .validationCheck("integrity.assignment", {
            .checkPartitionIntegrity(object)
            TRUE
        })
    } else {
        .validationUnknown("integrity.assignment",
            "Recomputation disabled or structural prerequisites failed")
    }
    out[[length(out) + 1L]] <- .validationUnknown("universe.verification",
        "Partition does not own the universe; validate its SegmentationSet")
    out[[length(out) + 1L]] <- .validationCheck("run_metadata.structure",
        if (methods::is(object@run_metadata, "SimpleList")) TRUE else
            "run_metadata must be a SimpleList")
    out[[length(out) + 1L]] <- .validationCheck("metadata.structure",
        if (methods::is(object@metadata, "SimpleList")) TRUE else
            "metadata must be a SimpleList")
    out[[length(out) + 1L]] <- .provenanceChecks(object)
    .validationBind(out)
}

#' Validate a segmentation partition explicitly
#'
#' Checks census identity, assignment conservation, ambiguity evidence and
#' assignment integrity independently of provenance completeness. Empty censuses
#' and zero-molecule cells are retained. External census completeness remains
#' unchecked; assignment mass cannot certify missing zero-molecule cells.
#' A partition does not own its universe;
#' standalone universe verification is UNKNOWN. No input is repaired or changed.
#' @param object A SegmentationPartition.
#' @param check_integrity Logical scalar. Recompute the existing assignment
#'   fingerprint by default. FALSE reports integrity as UNKNOWN. Fingerprints
#'   bind identity/assignments, not provenance, authenticity or spatial
#'   metadata.
#' @param strict Logical scalar. Stop only for FAIL/error findings, after
#'   checks.
#' @return A DataFrame with character check, status, severity and message
#'   columns.
#'   Statuses are PASS, FAIL, UNKNOWN and NOT_APPLICABLE; severities are error,
#'   warning and note. Unknown candidate evidence is reported without
#'   fabrication.
#' @examples
#' m <- S4Vectors::DataFrame(molecule_id = "m1", feature_id = "g1")
#' p <- SegmentationPartition(m, "s", S4Vectors::DataFrame(cell_id = "c1"),
#'     "m1", "assigned", "c1")
#' validateSegmentation(p, strict = TRUE)
#' @export
validateSegmentation <- function(object, check_integrity = TRUE,
                                 strict = FALSE) {
    .validationFlag(check_integrity)
    .validationFinish(.partitionValidation(object, check_integrity), strict)
}
