#' @rdname SegmentationPartition
#' @export
methods::setGeneric("cellCensus",
    function(x) methods::standardGeneric("cellCensus"))

#' @rdname SegmentationPartition
#' @export
methods::setMethod("cellCensus", "SegmentationPartition", function(x) x@cells)

#' @rdname SegmentationPartition
#' @export
methods::setGeneric("assignmentStatus",
    function(x) methods::standardGeneric("assignmentStatus"))

#' @rdname SegmentationPartition
#' @export
methods::setMethod("assignmentStatus", "SegmentationPartition", function(x) {
    c("assigned", "unassigned", "ambiguous")[x@status]
})

#' @rdname SegmentationPartition
#' @export
methods::setGeneric("ambiguityData",
    function(x) methods::standardGeneric("ambiguityData"))

#' @rdname SegmentationPartition
#' @export
methods::setMethod("ambiguityData", "SegmentationPartition",
    function(x) x@ambiguity)

#' @rdname SegmentationPartition
#' @export
methods::setGeneric("provenance",
    function(x) methods::standardGeneric("provenance"))

#' @rdname SegmentationPartition
#' @export
methods::setMethod("provenance", "SegmentationPartition",
    function(x) x@provenance)

#' @rdname SegmentationPartition
#' @export
methods::setGeneric("runMetadata",
    function(x) methods::standardGeneric("runMetadata"))

#' @rdname SegmentationPartition
#' @export
methods::setMethod("runMetadata", "SegmentationPartition",
    function(x) x@run_metadata)

#' @rdname SegmentationSet
#' @export
methods::setGeneric("segmentationNames",
    function(x) methods::standardGeneric("segmentationNames"))

#' @rdname SegmentationSet
#' @export
methods::setMethod("segmentationNames", "SegmentationSet",
    function(x) names(x@segmentations))

#' @rdname SegmentationSet
#' @export
methods::setGeneric("segmentationInfo",
    function(x) methods::standardGeneric("segmentationInfo"))

#' @rdname SegmentationSet
#' @export
methods::setMethod("segmentationInfo", "SegmentationSet", function(x) {
    p <- x@segmentations
    result <- S4Vectors::DataFrame(name = segmentationNames(x),
        n_cells = unname(vapply(p, function(z) nrow(z@cells), integer(1))),
        n_molecules = unname(vapply(p, function(z) length(z@status),
                                   integer(1))))
    # Direct column assignment avoids DataFrame's argument coercion machinery.
    result$run_metadata <- lapply(p, function(z) z@run_metadata)
    result
})

#' @rdname SegmentationSet
#' @export
methods::setGeneric("moleculeData",
    function(x) methods::standardGeneric("moleculeData"))

#' @rdname SegmentationSet
#' @export
methods::setMethod("moleculeData", "SegmentationSet", function(x) {
    m <- x@molecules
    result <- S4Vectors::DataFrame(molecule_id = m$molecule_id,
        feature_id = x@features$feature_id[m$feature_index],
        row.names = rownames(m))
    ancillary <- setdiff(colnames(m),
        c("molecule_id", "molecule_index", "feature_index"))
    for (column in ancillary) {
        result[[column]] <- m[[column]]
    }
    result
})

#' @rdname SegmentationSet
#' @export
methods::setGeneric("featureData",
    function(x) methods::standardGeneric("featureData"))

#' @rdname SegmentationSet
#' @export
methods::setMethod("featureData", "SegmentationSet", function(x) x@features)

#' @rdname SegmentationSet
#' @export
methods::setMethod("provenance", "SegmentationSet", function(x) x@provenance)
