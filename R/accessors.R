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
