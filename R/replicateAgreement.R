.replicateChoose2 <- function(n) {
    # Each count and their total are bounded by the integer-index universe.
    # Even n*(n-1) fits signed LLint at that limit; the sum fits choose(N,2).
    if (!.exactCounts(n) || sum(n) > .Machine$integer.max) {
        stop("Replicate counts exceed the supported integer-index range",
             call. = FALSE)
    }
    values <- methods::as(as.integer(n), "LLint")
    sum((values * (values - 1L)) %/% 2L)
}

.replicateARIFromCounts <- function(shared, row, column, total) {
    # Twice the adjusted denominator times T is
    # R*(T-C) + C*(T-R). Both products are nonnegative. Checking their
    # factors is exact without constructing products beyond signed LLint.
    if ((row == 0L || column == total) &&
        (column == 0L || row == total)) {
        return(NA_real_)
    }
    # Only the public ratio uses double. Exact differences retain small
    # discordant populations even when neighboring large counts round alike.
    # Audited decimal LLint conversion preserves the counts up to this final
    # ratio boundary; ordinary numeric parsing then deliberately rounds them.
    v <- as.numeric(as.character(c(shared, row - shared, column - shared,
        total - row - column + shared, row, column,
        total - column, total - row)))
    2 * (v[1L] * v[4L] - v[2L] * v[3L]) /
        (v[5L] * v[7L] + v[6L] * v[8L])
}

.replicateARI <- function(edges) {
    weight <- S4Vectors::mcols(edges)$weight
    n <- .countTotal(weight)
    if (is.na(n)) stop("Invalid replicate contingency counts", call. = FALSE)
    .checkIndexSize(n)
    if (n < 2L) return(NA_real_)
    row <- .sumByIndex(S4Vectors::from(edges), weight,
                       S4Vectors::nLnode(edges))
    column <- .sumByIndex(S4Vectors::to(edges), weight,
                          S4Vectors::nRnode(edges))
    .replicateARIFromCounts(.replicateChoose2(weight), .replicateChoose2(row),
        .replicateChoose2(column), .replicateChoose2(n))
}

.replicatePair <- function(p, r, name1, name2) {
    status <- .comparisonStatus(p, r)
    n <- as.numeric(length(p@status))
    agree <- sum(status$n_molecules[c(1L, 5L, 9L)])
    joint <- status$n_molecules[1L]
    edges <- .aggregatePairs(.sharedPairs(p, r), nrow(p@cells), nrow(r@cells))
    summary <- S4Vectors::DataFrame(segmentation_1 = name1,
        segmentation_2 = name2, n_molecules = n, n_status_agree = agree,
        n_status_disagree = n - agree,
        status_agreement_fraction = if (n > 0) agree / n else NA_real_,
        n_jointly_assigned = joint,
        jointly_assigned_fraction = if (n > 0) joint / n else NA_real_,
        ari_jointly_assigned = .replicateARI(edges))
    contingency <- S4Vectors::DataFrame(segmentation_1 = rep(name1, 9L),
        segmentation_2 = rep(name2, 9L), status)
    list(summary = summary, contingency = contingency)
}

#' Label-invariant agreement among segmentation runs
#'
#' Describe computational repeatability through full-universe assignment-status
#' accounting and conditional adjusted Rand index (ARI). These results do not
#' certify biological truth, segmentation accuracy or independent replication.
#' @md
#' @param object A valid SegmentationSet owning one immutable molecule universe.
#' @param segmentations A base character vector naming at least two distinct
#'   existing segmentations. Names must be nonmissing and nonempty. Supplied
#'   order determines pair order; there is no implicit selection or reference.
#' @details Every unordered selected-run pair is reported, with the earlier
#'   supplied run first: (1,2), (1,3), ..., (2,3), .... Construction has aligned
#'   assignments through lossless molecule identities. Analysis verifies the
#'   universe digest once and selected assignment digests once each, then uses
#'   those existing indices. Equal cell IDs in different runs have no matching
#'   semantics. Cell labels and census permutations do not change agreement.
#'
#'   Assignment-status agreement counts identical assigned, unassigned or
#'   ambiguous statuses over the complete universe. All nine status combinations
#'   are retained, including zeros. Ambiguous and unassigned molecules are never
#'   forced into cells or combined into a pseudo-cell. The jointly assigned
#'   count includes only molecules uniquely assigned in both runs. Both reported
#'   fractions use the full molecule universe, with NA for an empty universe.
#'
#'   ARI compares the partition co-membership of jointly assigned molecules,
#'   using a sparse shared-molecule contingency and its conditional margins.
#'   These margins do not replace any primary graph denominators. If A is the
#'   sum of choose(edge mass,2), R and C the corresponding row/column pair sums,
#'   and T=choose(jointly assigned count,2), the usual ARI is
#'   `(A-R*C/T)/((R+C)/2-R*C/T)`. Pair counts and degeneracy checks remain exact
#'   through signed LLint arithmetic at the supported integer-index universe
#'   limit. Zero denominator is detected through exact nonnegative factors,
#'   before double conversion for the final numeric ratio. ARI is NA_real_ for
#'   fewer than two jointly assigned molecules or an exactly zero adjusted
#'   denominator, including all-singleton and one-cell equivalent partitions.
#'
#'   All input censuses, assignment evidence, provenance and run metadata remain
#'   unchanged. Complete selected run metadata, including unknown or missing
#'   values, is retained without simplification. Integrity fingerprints do not
#'   establish provenance or authenticity. There is no consensus, averaged
#'   score, cell matching, geometry or biological-state integration. Graph,
#'   state, support and topology variation should be analyzed separately for
#'   each run.
#' @return A SimpleList containing exactly two DataFrames:
#'   * `pairwise_summary`: segmentation_1, segmentation_2, n_molecules,
#'     n_status_agree, n_status_disagree, status_agreement_fraction,
#'     n_jointly_assigned, jointly_assigned_fraction and ari_jointly_assigned.
#'   * `assignment_status_contingency`: segmentation_1, segmentation_2,
#'     status_1, status_2 and n_molecules. Each pair has nine rows in
#'     assigned, unassigned, ambiguous order, with side 1 varying slowest.
#'
#'   Names/statuses are character; molecule counts are exact whole-valued
#'   numeric values; fractions and ARI are doubles. Result metadata contains
#'   `segmentations` in supplied order and `run_metadata`, a named SimpleList
#'   retaining each selected partition's complete run metadata.
#' @seealso [SegmentationSet()], [compareSegmentations()], [assignmentSummary()]
#' @examples
#' m <- S4Vectors::DataFrame(molecule_id = paste0("m", 1:4),
#'                          feature_id = rep("g", 4))
#' a <- SegmentationPartition(m, "run1",
#'     S4Vectors::DataFrame(cell_id = c("x", "y", "empty")),
#'     m$molecule_id, rep("assigned", 4), c("x", "x", "y", "y"),
#'     run_metadata = S4Vectors::SimpleList(seed = 1L))
#' b <- SegmentationPartition(m, "run2",
#'     S4Vectors::DataFrame(cell_id = c("unrelated", "labels")),
#'     m$molecule_id, rep("assigned", 4),
#'     c("unrelated", "unrelated", "labels", "labels"),
#'     run_metadata = S4Vectors::SimpleList(seed = 2L))
#' ss <- SegmentationSet(m, list(run1 = a, run2 = b), "specimen")
#' result <- replicateAgreement(ss, c("run1", "run2"))
#' result$pairwise_summary
#' result$assignment_status_contingency
#' S4Vectors::metadata(result)$run_metadata
#' @export
replicateAgreement <- function(object, segmentations) {
    if (!methods::is(object, "SegmentationSet")) {
        stop("object must be a SegmentationSet", call. = FALSE)
    }
    if (!identical(class(segmentations), "character") ||
        !is.null(dim(segmentations)) || length(segmentations) < 2L ||
        anyNA(segmentations) || any(!nzchar(segmentations)) ||
        anyDuplicated(segmentations) ||
        !all(segmentations %in% names(object@segmentations))) {
        stop("segmentations must name at least two distinct existing runs",
             call. = FALSE)
    }
    methods::validObject(object)
    .validateUniverse(S4Vectors::SimpleList(molecules = object@molecules,
        features = object@features,
        fingerprint_value = object@universe_fingerprint,
        fingerprint_version = object@fingerprint_version,
        fingerprint_algorithm = object@fingerprint_algorithm))
    runs <- object@segmentations[segmentations]
    for (i in seq_along(runs)) .checkPartitionIntegrity(runs[[i]])
    nr <- length(runs)
    np <- as.double(nr) * (nr - 1L) / 2
    .checkIndexSize(9 * np)
    pairs <- vector("list", np)
    k <- 0L
    for (i in seq_len(nr - 1L)) {
        for (j in seq.int(i + 1L, nr)) {
            k <- k + 1L
            pairs[[k]] <- .replicatePair(runs[[i]], runs[[j]],
                segmentations[[i]], segmentations[[j]])
        }
    }
    result <- S4Vectors::SimpleList(
        pairwise_summary = do.call(rbind, lapply(pairs, `[[`, "summary")),
        assignment_status_contingency = do.call(rbind,
            lapply(pairs, `[[`, "contingency")))
    source_runs <- S4Vectors::SimpleList(lapply(seq_along(runs),
        function(i) runs[[i]]@run_metadata))
    names(source_runs) <- segmentations
    S4Vectors::metadata(result) <- list(segmentations = segmentations,
                                      run_metadata = source_runs)
    result
}
