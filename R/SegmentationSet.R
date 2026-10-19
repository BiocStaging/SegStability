# Assisted-by: OpenAI Codex
.partitionBindingProblems <- function(partition, n, fingerprint, version,
                                      algorithm) {
    if (!methods::is(partition, "SegmentationPartition")) {
        return("Every segmentation must be a SegmentationPartition")
    }
    problems <- methods::validObject(partition, test = TRUE)
    if (!isTRUE(problems)) {
        return(problems)
    }
    if (length(partition@status) != n ||
        !identical(partition@universe_fingerprint, fingerprint) ||
        !identical(partition@fingerprint_version, version) ||
        !identical(partition@fingerprint_algorithm, algorithm)) {
        return("Partition universe binding mismatch")
    }
    TRUE
}

.setValidity <- function(object) {
    if (!.scalarString(object@specimen_id)) {
        return("specimen_id must be one nonmissing nonempty character string")
    }
    universe <- tryCatch({
        .checkUniverseStructure(object@molecules, object@features)
        TRUE
    }, error = function(e) conditionMessage(e))
    if (!isTRUE(universe)) {
        return(universe)
    }
    if (!.scalarString(object@universe_fingerprint) ||
        !grepl("^[0-9a-f]{64}$", object@universe_fingerprint) ||
        !identical(object@fingerprint_version, "1") ||
        !identical(object@fingerprint_algorithm, "sha256") ||
        !identical(object@schema_version, "1")) {
        return("Invalid set fingerprint metadata or schema_version")
    }
    partitions <- object@segmentations
    nm <- names(partitions)
    if (!length(partitions) || is.null(nm) || anyNA(nm) ||
        any(!nzchar(nm)) || anyDuplicated(nm)) {
        return("segmentations require at least one uniquely named partition")
    }
    for (i in seq_along(partitions)) {
        p <- partitions[[i]]
        binding <- .partitionBindingProblems(p, nrow(object@molecules),
            object@universe_fingerprint, object@fingerprint_version,
            object@fingerprint_algorithm)
        if (!isTRUE(binding)) {
            return(binding)
        }
        if (!identical(nm[i], p@name)) {
            return("Collection names must equal partition names")
        }
    }
    TRUE
}

#' A named collection sharing one immutable molecule universe
#'
#' Own the molecule universe and feature dictionary once and bind independently
#' censused segmentations to it. This is descriptive structural infrastructure,
#' not a segmentation quality score or ranking.
#' @rdname SegmentationSet
#' @md
#' @details Construction preserves molecule order and lossless identity using
#'   the certified ID paths documented in [SegmentationPartition()]. Feature
#'   IDs are base character; the dictionary follows first occurrence. Partitions
#'   contain aligned assignment indices and binding fingerprints, without copies
#'   of molecule IDs or the feature dictionary. Full independent cell censuses,
#'   including zero-molecule cells, remain unchanged.
#'
#'   At least one partition is required. Collection names must be unique,
#'   nonmissing, nonempty and identical to partition names. Input order is
#'   preserved. Empty molecule universes and empty cell censuses are supported.
#'
#'   The universe is prepared and fingerprinted once. Each partition's stored
#'   molecule count, universe fingerprint, version and algorithm must match;
#'   its assignment fingerprint is checked for integrity. Universe order, ID
#'   family and feature identity/mapping are integrity-significant. Version 1
#'   SHA-256 fingerprints do not establish provenance or source authenticity.
#'   Ancillary metadata is outside the identity/assignment digest scope.
#'
#'   S4 validity checks structure and binding metadata without recomputing
#'   digests. The universe is scientifically immutable; do not modify core
#'   slots.
#'   There are no public replacement methods. Ordinary R copy-on-modify applies
#'   to returned DataFrames and SimpleLists. Reference-valued ancillary metadata
#'   retains its normal R semantics; the package does not deep-copy
#'   environments.
#'
#'   Coordinate information records declarations only. Construction does not
#'   certify spatial compatibility or transform coordinates. Unknown provenance
#'   is retained explicitly. Serialization uses standard R objects; no cross-R
#'   binary compatibility guarantee is added.
#' @seealso [SegmentationPartition()]
#' @param molecules A data.frame or DataFrame with unique lossless molecule_id
#'   and character feature_id, in immutable universe order. Other columns are
#'   retained.
#' @param segmentations A uniquely named base list or S4Vectors SimpleList of
#'   valid SegmentationPartition objects, with matching names and universe
#'   binding.
#' @param specimen_id Exactly one nonmissing, nonempty character string.
#' @param coordinate_info S4Vectors SimpleList of coordinate declarations.
#' @param provenance S4Vectors SimpleList of provenance. Omitted standardized
#'   concepts receive explicit unknown placeholders, as for partitions.
#' @param metadata S4Vectors SimpleList of ancillary set metadata.
#' @param x A SegmentationSet.
#' @param object A SegmentationSet to display or summarize.
#' @param segmentation One exact existing segmentation name.
#' @param i One exact nonempty segmentation name for retrieval.
#' @param j Not supported.
#' @param ... Additional arguments to summary (unused), or only
#'   \code{exact = TRUE} for retrieval.
#' @return SegmentationSet() returns a validated SegmentationSet.
#'   segmentationNames() returns names in collection order. segmentationInfo()
#'   returns a DataFrame with name, n_cells, n_molecules and a run_metadata list
#'   column retaining each partition's complete SimpleList without
#'   simplification.
#'   moleculeData() returns a DataFrame with external molecule_id, reconstructed
#'   feature_id and ancillary columns, without internal indices. featureData()
#'   returns the feature dictionary DataFrame. provenance() returns set
#'   provenance.
#'   \code{x[[name]]} returns the existing universe-bound partition after
#'   integrity checking. Only one exact character name is supported; partial
#'   or positional selection and replacement are not supported.
#'   assignmentData(x, segmentation) returns exactly molecule_id, semantic
#'   status and native typed cell_id in owning-universe order, one row per
#'   molecule. Unassigned and ambiguous cell IDs are typed missing; candidate
#'   evidence remains in ambiguityData(). Integrity is checked on read.
#'   summary() returns a SimpleList with specimen_id, n_molecules, n_features,
#'   n_segmentations, segmentations (the segmentationInfo() DataFrame), and
#'   coordinate_info (the exact stored declaration, including unknowns, without
#'   validation, completion or spatial prerequisites).
#'   show() displays structural counts and names. Neither method compares cells.
#'   runMetadata(), cellCensus(), assignmentStatus() and ambiguityData() remain
#'   partition accessors.
#' @slot molecules Immutable molecule DataFrame with internal positional
#'   indices.
#' @slot features First-occurrence feature dictionary DataFrame.
#' @slot segmentations Named SimpleList of universe-bound partitions.
#' @slot specimen_id Specimen identity.
#' @slot coordinate_info Coordinate declarations SimpleList.
#' @slot provenance Set provenance SimpleList.
#' @slot universe_fingerprint Shared universe identity digest.
#' @slot fingerprint_version Canonical fingerprint serialization version.
#' @slot fingerprint_algorithm Digest algorithm.
#' @slot schema_version Object schema version.
#' @slot metadata Ancillary metadata SimpleList.
#' @examples
#' molecules <- S4Vectors::DataFrame(
#'     molecule_id = c("m1", "m2"), feature_id = c("g", "h"))
#' p <- SegmentationPartition(molecules, "native",
#'     S4Vectors::DataFrame(cell_id = c("c1", "unused")),
#'     molecule_id = c("m2", "m1"),
#'     status = c("unassigned", "assigned"),
#'     cell_id = c(NA_character_, "c1"))
#' ss <- SegmentationSet(molecules, list(native = p), "specimen1")
#' ss
#' segmentationNames(ss)
#' ss[["native"]]
#' assignmentData(ss, "native")
#' segmentationInfo(ss)
#' moleculeData(ss)
#' featureData(ss)
#' provenance(ss)
#' summary(ss)
#' @exportClass SegmentationSet
methods::setClass("SegmentationSet", slots = c(
    molecules = "DataFrame", features = "DataFrame",
    segmentations = "SimpleList",
    specimen_id = "character", coordinate_info = "SimpleList",
    provenance = "SimpleList", universe_fingerprint = "character",
    fingerprint_version = "character", fingerprint_algorithm = "character",
    schema_version = "character", metadata = "SimpleList"),
    validity = .setValidity)

#' @rdname SegmentationSet
#' @export
SegmentationSet <- function(molecules, segmentations, specimen_id,
    coordinate_info = S4Vectors::SimpleList(),
    provenance = S4Vectors::SimpleList(), metadata = S4Vectors::SimpleList()) {
    if (!.scalarString(specimen_id)) {
        stop("specimen_id must be one nonmissing nonempty character string",
             call. = FALSE)
    }
    if (!methods::is(segmentations, "SimpleList")) {
        if (!identical(class(segmentations), "list")) {
            stop("segmentations must be a named list or SimpleList",
                 call. = FALSE)
        }
        segmentations <- S4Vectors::SimpleList(segmentations)
    }
    if (!methods::is(coordinate_info, "SimpleList") ||
        !methods::is(metadata, "SimpleList")) {
        stop("coordinate_info and metadata must be SimpleList", call. = FALSE)
    }
    provenance <- .completeProvenance(provenance)
    universe <- .prepareUniverse(molecules)
    object <- methods::new("SegmentationSet", molecules = universe$molecules,
        features = universe$features, segmentations = segmentations,
        specimen_id = specimen_id, coordinate_info = coordinate_info,
        provenance = provenance,
        universe_fingerprint = universe$fingerprint_value,
        fingerprint_version = universe$fingerprint_version,
        fingerprint_algorithm = universe$fingerprint_algorithm,
        schema_version = "1", metadata = metadata)
    for (i in seq_along(object@segmentations)) {
        .checkPartitionIntegrity(object@segmentations[[i]])
    }
    object
}

.checkSetIntegrity <- function(object) {
    methods::validObject(object)
    .validateUniverse(S4Vectors::SimpleList(molecules = object@molecules,
        features = object@features,
        fingerprint_value = object@universe_fingerprint,
        fingerprint_version = object@fingerprint_version,
        fingerprint_algorithm = object@fingerprint_algorithm))
    for (i in seq_along(object@segmentations)) {
        .checkPartitionIntegrity(object@segmentations[[i]])
    }
    invisible(NULL)
}

#' @rdname SegmentationSet
#' @export
methods::setMethod("show", "SegmentationSet", function(object) {
    cat("SegmentationSet:", object@specimen_id, "\n")
    cat(" ", nrow(object@molecules), "molecules;", nrow(object@features),
        "features;", length(object@segmentations), "segmentations\n")
    cat(" ", paste(segmentationNames(object), collapse = ", "), "\n")
})

#' @rdname SegmentationSet
#' @export
methods::setMethod("summary", "SegmentationSet", function(object, ...) {
    S4Vectors::SimpleList(specimen_id = object@specimen_id,
        n_molecules = nrow(object@molecules),
        n_features = nrow(object@features),
        n_segmentations = length(object@segmentations),
        segmentations = segmentationInfo(object),
        coordinate_info = object@coordinate_info)
})
