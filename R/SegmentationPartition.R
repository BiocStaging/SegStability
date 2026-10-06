.scalarString <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.ambiguityProblems <- function(ambiguity, status, ncell) {
    if (!identical(sort(names(ambiguity)), c("candidates", "evidence"))) {
        return("ambiguity must contain evidence and candidates")
    }
    e <- ambiguity$evidence
    h <- ambiguity$candidates
    if (!methods::is(e, "DataFrame") || !methods::is(h, "Hits") ||
        !all(c("molecule_index", "candidate_count", "candidate_ids_known") %in%
             colnames(e))) {
        return("Invalid ambiguity evidence or candidates structure")
    }
    if (!identical(e$molecule_index, which(status == 3L)) ||
        !is.integer(e$candidate_count) ||
        any(e$candidate_count < 0L, na.rm = TRUE) ||
        !is.logical(e$candidate_ids_known) || anyNA(e$candidate_ids_known)) {
        return("Invalid ambiguity indices, candidate_count or known flags")
    }
    if (S4Vectors::nLnode(h) != length(status) ||
        S4Vectors::nRnode(h) != ncell) {
        return("Ambiguity candidate dimensions disagree with universe/census")
    }
    from <- S4Vectors::from(h)
    to <- S4Vectors::to(h)
    if (anyNA(from) || anyNA(to) || any(from < 1L | from > length(status)) ||
        any(to < 1L | to > ncell) || any(!from %in% e$molecule_index) ||
        anyDuplicated(paste(from, to, sep = ":")) ||
        !identical(base::order(from, to, method = "radix"), seq_along(from))) {
        return("Invalid or noncanonical ambiguity candidate references")
    }
    observed <- tabulate(from, nbins = length(status))[e$molecule_index]
    known <- e$candidate_ids_known
    if (any(observed[!known] != 0L) || anyNA(e$candidate_count[known]) ||
        any(e$candidate_count[known] != observed[known])) {
        return("Ambiguity candidate counts/known flags are inconsistent")
    }
    TRUE
}

.partitionValidity <- function(object) {
    if (!.scalarString(object@name)) {
        return("name must be a nonempty nonmissing character scalar")
    }
    census <- tryCatch({
        .validateIds(object@cells$cell_id, "cell census")
        TRUE
    }, error = function(e) conditionMessage(e))
    if (!isTRUE(census)) {
        return(census)
    }
    s <- object@status
    i <- object@cell_index
    if (length(s) != length(i) || anyNA(s) || any(!s %in% 1:3)) {
        return("status/cell_index lengths or status codes are invalid")
    }
    assigned <- s == 1L
    if (anyNA(i[assigned]) ||
        any(i[assigned] < 1L | i[assigned] > nrow(object@cells)) ||
        any(!is.na(i[!assigned]))) {
        return("Invalid assigned cell_index or forced unresolved assignment")
    }
    ambiguity <- .ambiguityProblems(object@ambiguity, s, nrow(object@cells))
    if (!isTRUE(ambiguity)) {
        return(ambiguity)
    }
    hashes <- c(object@universe_fingerprint, object@assignment_fingerprint)
    if (length(object@universe_fingerprint) != 1L ||
        length(object@assignment_fingerprint) != 1L || anyNA(hashes) ||
        any(!grepl("^[0-9a-f]{64}$", hashes)) ||
        !identical(object@fingerprint_version, "1") ||
        !identical(object@fingerprint_algorithm, "sha256") ||
        !identical(object@schema_version, "1")) {
        return(paste("Invalid fingerprint fields, version/algorithm",
                     "or schema_version"))
    }
    TRUE
}

#' One segmentation bound to an immutable molecule universe
#'
#' Represent a segmentation with an explicit cell census, exactly aligned
#' molecule assignments, and retained ambiguity evidence. These objects support
#' descriptive sensitivity analysis; they do not score segmentation accuracy.
#' @rdname SegmentationPartition
#'
#' @details Molecule and cell identifiers support base integer and character,
#'   and certified exact paths for S4Vectors LLint and optional bit64 integer64.
#'   Ordinary double, logical, factor, complex, raw, missing and empty
#'   identifiers
#'   are rejected. Matching requires the same family and never converts IDs
#'   through ordinary double. Original cell IDs and census order are preserved,
#'   including cells with no molecules. No cell ID is an implicit sentinel.
#'
#'   For v0.99, feature IDs support base character only. Missing and empty
#'   values
#'   fail; whitespace, spelling, case and Unicode sequences are preserved
#'   without
#'   normalization. There is no automatic conversion from other types. Repeated
#'   features across molecules are valid; internal dictionary entries are unique
#'   in first-occurrence order and map deterministically to integer feature
#'   indices.
#'
#'   Construction verifies exact molecule membership and reorders assignments
#'   once to the supplied universe order. Assignment-side molecule IDs are then
#'   discarded. Internal status codes and positional indices are implementation
#'   details; use the accessors. Do not mutate scientific core slots.
#'
#'   Ambiguity input uses external IDs. Its evidence table has molecule_id,
#'   candidate_count (nonnegative whole count or NA), and candidate_ids_known
#'   (nonmissing logical). Its candidates table has molecule_id and cell_id.
#'   Both tables use the same ID families as their universe/census. Candidate
#'   pairs must be unique. Known sets are complete and counts must agree;
#'   omitted known counts are derived from the supplied set. Count-only evidence
#'   has a known count and FALSE known flag; entirely unknown evidence has NA
#'   count and FALSE flag. Omitted evidence for ambiguous molecules remains
#'   entirely unknown. No candidate identities are fabricated. Counts do not
#'   impose an additional biological threshold on unresolved assignments.
#'
#'   Internal version 1 SHA-256 fingerprints bind ordered molecule/feature
#'   identity
#'   and ordered census/assignment/ambiguity content. Canonical encoding uses
#'   framed fields, exact decimal integer keys, ID-family tags and missingness
#'   tags. Character text is encoded in UTF-8 without normalization; byte-marked
#'   strings retain their bytes. Census and universe order affect the digest;
#'   candidate input order does not. Optional molecule/cell metadata, run
#'   metadata
#'   and provenance are outside this identity/assignment digest scope.
#'   Fingerprints
#'   establish neither source authenticity nor provenance. Structural validity
#'   checks fingerprint fields without repeating expensive source-ID matching.
#'
#' @param molecules A S4Vectors DataFrame with unique, lossless molecule_id and
#'   base character feature_id columns. Row order defines the immutable
#'   universe.
#'   Additional columns are retained internally without transformation. Internal
#'   molecule_index and feature_index columns are reserved.
#' @param name Nonempty, nonmissing character scalar naming the segmentation.
#' @param cells A S4Vectors DataFrame containing unique, lossless cell_id values
#'   for the full census, including zero-molecule cells.
#' @param molecule_id External assignment molecule IDs, exactly the universe
#'   members in any order, using the universe ID family.
#' @param status Character vector parallel to molecule_id, containing only
#'   assigned, unassigned or ambiguous.
#' @param cell_id Cell IDs parallel to molecule_id. Assigned entries must name
#'   census cells; unresolved entries must be missing in the census ID family.
#'   NULL is permitted only when there are no assigned molecules.
#' @param ambiguity NULL for unknown candidate evidence, or a S4Vectors
#'   SimpleList
#'   containing external-ID-keyed evidence and candidates DataFrames.
#' @param run_metadata S4Vectors SimpleList of run metadata; empty by default.
#' @param provenance S4Vectors SimpleList of provenance. Missing standardized
#'   concepts receive explicit unknown placeholders; supplied values are
#'   retained.
#' @param metadata S4Vectors SimpleList of ancillary metadata; empty by default.
#' @param x A SegmentationPartition.
#' @param object A SegmentationPartition to display.
#' @return SegmentationPartition() returns a validated SegmentationPartition.
#'   cellCensus() returns the full cell DataFrame; assignmentStatus() returns
#'   semantic character statuses in universe order; ambiguityData() returns a
#'   SimpleList with index-based evidence and molecule-to-cell Hits.
#'   provenance() and runMetadata() return their respective SimpleLists.
#' @slot name Segmentation name.
#' @slot cells Explicit lossless cell census DataFrame.
#' @slot status Internal integer status enum in universe order.
#' @slot cell_index Integer census row for assigned molecules, NA otherwise.
#' @slot ambiguity SimpleList containing evidence DataFrame and candidate Hits.
#' @slot run_metadata Structured run metadata SimpleList.
#' @slot provenance Structured provenance SimpleList, retaining unknowns.
#' @slot universe_fingerprint Universe identity digest as character.
#' @slot assignment_fingerprint Assignment identity digest as character.
#' @slot fingerprint_version Canonical serialization version as character.
#' @slot fingerprint_algorithm Digest algorithm as character.
#' @slot schema_version Object schema version as character.
#' @slot metadata Ancillary metadata SimpleList.
#' @examples
#' molecules <- S4Vectors::DataFrame(
#'     molecule_id = c("m1", "m2", "m3"), feature_id = c("g", "g", "h"))
#' cells <- S4Vectors::DataFrame(cell_id = c("c1", "unused"))
#' p <- SegmentationPartition(molecules, "native", cells,
#'     molecule_id = c("m3", "m1", "m2"),
#'     status = c("ambiguous", "assigned", "unassigned"),
#'     cell_id = c(NA_character_, "c1", NA_character_))
#' p
#' cellCensus(p)
#' assignmentStatus(p)
#' ambiguityData(p)
#' provenance(p)
#' runMetadata(p)
#' @importClassesFrom S4Vectors DataFrame SimpleList Hits
#' @exportClass SegmentationPartition
methods::setClass("SegmentationPartition", slots = c(
    name = "character", cells = "DataFrame", status = "integer",
    cell_index = "integer", ambiguity = "SimpleList",
    run_metadata = "SimpleList", provenance = "SimpleList",
    universe_fingerprint = "character", assignment_fingerprint = "character",
    fingerprint_version = "character", fingerprint_algorithm = "character",
    schema_version = "character", metadata = "SimpleList"),
    validity = .partitionValidity)

.completeProvenance <- function(provenance) {
    if (!methods::is(provenance, "SimpleList") ||
        (length(provenance) && (is.null(names(provenance)) ||
         anyNA(names(provenance)) || any(!nzchar(names(provenance))) ||
         anyDuplicated(names(provenance))))) {
        stop("provenance must be a uniquely named SimpleList", call. = FALSE)
    }
    fields <- c("source", "source_version", "retrieval_date", "source_hashes",
        "software", "software_version", "software_hashes", "command",
        "parameters", "model", "threads", "runtime", "observation_window",
        "output_hashes")
    for (field in setdiff(fields, names(provenance))) {
        provenance[[field]] <- switch(field,
            threads = NA_integer_, runtime = NA_real_,
            parameters = S4Vectors::SimpleList(),
            observation_window = S4Vectors::SimpleList(), NA_character_)
    }
    provenance
}

#' @rdname SegmentationPartition
#' @export
SegmentationPartition <- function(molecules, name, cells, molecule_id, status,
    cell_id = NULL, ambiguity = NULL, run_metadata = S4Vectors::SimpleList(),
    provenance = S4Vectors::SimpleList(), metadata = S4Vectors::SimpleList()) {
    if (!.scalarString(name)) {
        stop("name must be a nonempty nonmissing character scalar",
             call. = FALSE)
    }
    universe <- .prepareUniverse(molecules)
    if (!methods::is(cells, "DataFrame") ||
        !"cell_id" %in% colnames(cells) || anyDuplicated(colnames(cells))) {
        stop("cells must be a DataFrame containing cell_id", call. = FALSE)
    }
    .checkIndexSize(nrow(cells))
    .validateIds(cells$cell_id, "cell census")
    order <- .alignMolecules(molecule_id, universe$molecules$molecule_id)
    status <- .canonicalStatus(status)
    if (length(status) != length(molecule_id) ||
        (!is.null(cell_id) && length(cell_id) != length(molecule_id))) {
        stop("Assignment vector lengths must agree", call. = FALSE)
    }
    status <- status[order]
    if (!is.null(cell_id)) {
        cell_id <- cell_id[order]
    }
    cell_index <- .mapCells(cell_id, cells$cell_id, status)
    ambiguity <- .prepareAmbiguity(ambiguity, universe$molecules$molecule_id,
                                   cells$cell_id, status)
    provenance <- .completeProvenance(provenance)
    if (!methods::is(run_metadata, "SimpleList") ||
        !methods::is(metadata, "SimpleList")) {
        stop("run_metadata and metadata must be SimpleList", call. = FALSE)
    }
    fingerprint <- .assignmentFingerprint(universe$fingerprint_value, cells,
                                           status, cell_index, ambiguity)
    object <- methods::new("SegmentationPartition", name = name, cells = cells,
        status = status, cell_index = cell_index, ambiguity = ambiguity,
        run_metadata = run_metadata, provenance = provenance,
        universe_fingerprint = universe$fingerprint_value,
        assignment_fingerprint = fingerprint$fingerprint_value,
        fingerprint_version = fingerprint$fingerprint_version,
        fingerprint_algorithm = fingerprint$fingerprint_algorithm,
        schema_version = "1", metadata = metadata)
    methods::validObject(object)
    object
}

.checkPartitionIntegrity <- function(object, universe = NULL) {
    methods::validObject(object)
    if (!is.null(universe)) {
        .validateUniverse(universe)
        if (length(object@status) != nrow(universe$molecules) ||
            !identical(object@universe_fingerprint,
                       universe$fingerprint_value) ||
            !identical(object@fingerprint_version,
                       universe$fingerprint_version) ||
            !identical(object@fingerprint_algorithm,
                       universe$fingerprint_algorithm)) {
            stop("Partition universe binding mismatch", call. = FALSE)
        }
    }
    calculated <- .assignmentFingerprint(object@universe_fingerprint,
        object@cells, object@status, object@cell_index, object@ambiguity)
    if (!identical(object@assignment_fingerprint,
                   calculated$fingerprint_value)) {
        stop("Partition assignment fingerprint mismatch", call. = FALSE)
    }
    invisible(NULL)
}

#' @rdname SegmentationPartition
#' @export
methods::setMethod("show", "SegmentationPartition", function(object) {
    counts <- tabulate(object@status, nbins = 3L)
    cat("SegmentationPartition:", object@name, "\n")
    cat(" ", nrow(object@cells), "cells;", length(object@status), "molecules\n")
    cat("  assigned:", counts[1L], "unassigned:", counts[2L],
        "ambiguous:", counts[3L], "\n")
})
