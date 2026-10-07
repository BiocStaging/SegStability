.sharedPairs <- function(p1, p2) {
    # Alignment was certified at construction; these are positional indices.
    both <- which(p1@status == 1L & p2@status == 1L)
    S4Vectors::DataFrame(molecule_index = both, query = p1@cell_index[both],
        subject = p2@cell_index[both], weight = rep(1, length(both)))
}

.aggregatePairs <- function(pairs, n1, n2) {
    q <- pairs$query
    s <- pairs$subject
    o <- order(q, s, method = "radix")
    q <- q[o]
    s <- s[o]
    if (length(q)) {
        start <- which(c(TRUE, q[-1L] != q[-length(q)] |
                                  s[-1L] != s[-length(s)]))
        weight <- as.numeric(diff(c(start, length(q) + 1)))
        q <- q[start]
        s <- s[start]
    } else {
        weight <- numeric()
    }
    if (!.exactCounts(weight)) stop("Shared counts exceed exact count range")
    h <- S4Vectors::Hits(q, s, n1, n2)
    S4Vectors::mcols(h) <- S4Vectors::DataFrame(weight = weight)
    h
}

.directionalDominance <- function(index, weight, n) {
    maximum <- numeric(n)
    for (k in seq_along(index)) {
        i <- index[k]
        if (weight[k] > maximum[i]) maximum[i] <- weight[k]
    }
    # Dominance maximizes w/N. N is constant across each endpoint's edges,
    # so exact weight comparisons preserve fraction maxima and every exact tie.
    dominant <- weight == maximum[index]
    count <- tabulate(index[dominant], n)
    list(dominant = dominant, unique = dominant & count[index] == 1L,
         count = count)
}

.graphComponents <- function(query, subject, n1, n2) {
    n <- as.double(n1) + n2
    if (n > .Machine$integer.max) {
        stop("Combined census exceeds supported positional index range",
             call. = FALSE)
    }
    parent <- seq_len(n)
    rank <- integer(n)
    root <- function(i) {
        r <- i
        while (parent[r] != r) r <- parent[r]
        while (parent[i] != i) {
            next_i <- parent[i]
            parent[i] <<- r
            i <- next_i
        }
        r
    }
    for (k in seq_along(query)) {
        a <- root(query[k])
        b <- root(n1 + subject[k])
        if (a == b) next
        if (rank[a] < rank[b]) {
            parent[a] <- b
        } else {
            parent[b] <- a
            if (rank[a] == rank[b]) rank[a] <- rank[a] + 1L
        }
    }
    roots <- vapply(seq_len(n), root, integer(1))
    # First occurrence fixes IDs independently of union root choices.
    membership <- match(roots, unique(roots))
    c1 <- membership[seq_len(n1)]
    c2 <- membership[n1 + seq_len(n2)]
    nc <- length(unique(roots))
    sizes1 <- tabulate(c1, nc)
    sizes2 <- tabulate(c2, nc)
    list(side1 = c1, side2 = c2,
        table = S4Vectors::DataFrame(component_id = seq_len(nc),
            n_segmentation_1 = sizes1, n_segmentation_2 = sizes2,
            component_class = .componentClass(sizes1, sizes2)))
}

.assignedMass <- function(partition) {
    as.numeric(tabulate(partition@cell_index[partition@status == 1L],
                        nrow(partition@cells)))
}

.comparisonCells <- function(partition, assigned, index, weight, dominance,
                             component) {
    n <- nrow(partition@cells)
    matched <- .sumByIndex(index, weight, n)
    fraction <- rep(NA_real_, n)
    positive <- assigned > 0
    fraction[positive] <- matched[positive] / assigned[positive]
    S4Vectors::DataFrame(cell_id = partition@cells$cell_id,
        assigned_mass = assigned, matched_mass = matched,
        unmatched_mass = assigned - matched,
        partner_count = tabulate(index, n), dominance_count = dominance$count,
        component_id = component, matched_fraction = fraction,
        row.names = rownames(partition@cells))
}

.comparisonStatus <- function(p1, p2) {
    statuses <- c("assigned", "unassigned", "ambiguous")
    S4Vectors::DataFrame(status_1 = rep(statuses, each = 3L),
        status_2 = rep(statuses, 3L),
        n_molecules = as.numeric(tabulate((p1@status - 1L) * 3L + p2@status,
                                         9L)))
}

.compareAlignedPartitions <- function(object, nm, p, r) {
    # Independent denominators precede shared-pair construction and selection.
    n1 <- .assignedMass(p)
    n2 <- .assignedMass(r)
    h <- .aggregatePairs(.sharedPairs(p, r), nrow(p@cells), nrow(r@cells))
    q <- S4Vectors::from(h)
    s <- S4Vectors::to(h)
    w <- S4Vectors::mcols(h)$weight
    d1 <- .directionalDominance(q, w, nrow(p@cells))
    d2 <- .directionalDominance(s, w, nrow(r@cells))
    comp <- .graphComponents(q, s, nrow(p@cells), nrow(r@cells))
    c1 <- .comparisonCells(p, n1, q, w, d1, comp$side1)
    c2 <- .comparisonCells(r, n2, s, w, d2, comp$side2)
    S4Vectors::mcols(h) <- S4Vectors::DataFrame(weight = w,
        fraction_1 = w / c1$assigned_mass[q],
        fraction_2 = w / c2$assigned_mass[s],
        dominant_1 = d1$dominant, dominant_2 = d2$dominant,
        unique_dominant_1 = d1$unique, unique_dominant_2 = d2$unique,
        reciprocal_dominant = d1$unique & d2$unique)
    source_provenance <- S4Vectors::SimpleList(
        list(p@provenance, r@provenance))
    source_runs <- S4Vectors::SimpleList(list(p@run_metadata, r@run_metadata))
    names(source_provenance) <- names(source_runs) <- nm
    fingerprints <- c(p@assignment_fingerprint, r@assignment_fingerprint)
    names(fingerprints) <- nm
    .newComparison(segmentation_names = nm, cell_summary_1 = c1,
        cell_summary_2 = c2, edges = h, components = comp$table,
        assignment_summary = .comparisonStatus(p, r),
        universe_size = as.numeric(nrow(object@molecules)),
        universe_fingerprint = object@universe_fingerprint,
        assignment_fingerprints = fingerprints,
        fingerprint_version = object@fingerprint_version,
        fingerprint_algorithm = object@fingerprint_algorithm,
        schema_version = "1", coordinate_info = object@coordinate_info,
        provenance = .completeProvenance(S4Vectors::SimpleList(
            set = object@provenance, segmentations = source_provenance)),
        metadata = S4Vectors::SimpleList(specimen_id = object@specimen_id,
                                       run_metadata = source_runs))
}

#' Compare two segmentations of one immutable molecule universe
#'
#' Construct the complete positive shared-molecule correspondence graph and
#' descriptive structural summaries before any selection.
#' @md
#' @param object A valid SegmentationSet owning the shared molecule universe.
#' @param segmentation1,segmentation2 Distinct existing segmentation names,
#'   each one nonmissing, nonempty character string. Order determines query
#'   and subject sides and split/merge direction.
#' @details Both source assignments must retain valid shared-universe bindings
#'   and assignment fingerprints. The universe digest is checked once and both
#'   selected assignment digests are checked without repeating source-ID
#'   matching. Only uniquely assigned molecules on both sides contribute to
#'   edges; each supplies one unit to one pair before aggregation. Each own-side
#'   denominator counts all its uniquely assigned molecules, including those
#'   unresolved on the opposite side. Full censuses and isolates are retained.
#'
#'   Empty censuses/universes and empty positive graphs are supported. Counts
#'   must be exact whole values at most 2^53; scientific identifiers are never
#'   coerced to ordinary double. No cell-by-cell dense matrix is constructed.
#'   Combined census sizes beyond the integer positional-index range fail.
#'   No support, state, geometry or anchor selection is performed.
#' @return A validated SegmentationComparison. See
#'   [SegmentationComparison-class]
#'   for graph columns, side summaries, status accounting and component classes.
#' @seealso [correspondenceEdges()], [cellSummary()], [componentSummary()],
#'   [assignmentSummary()]
#' @examples
#' m <- S4Vectors::DataFrame(molecule_id = "m", feature_id = "g")
#' p <- SegmentationPartition(m, "first",
#'     S4Vectors::DataFrame(cell_id = "x"), "m", "assigned", "x")
#' q <- SegmentationPartition(m, "second",
#'     S4Vectors::DataFrame(cell_id = "a"), "m", "unassigned", NA_character_)
#' ss <- SegmentationSet(m, list(first = p, second = q), "specimen")
#' cmp <- compareSegmentations(ss, "first", "second")
#' cellSummary(cmp, "first")
#' assignmentSummary(cmp)
#' @export
compareSegmentations <- function(object, segmentation1, segmentation2) {
    if (!methods::is(object, "SegmentationSet")) {
        stop("object must be a SegmentationSet", call. = FALSE)
    }
    if (!.scalarString(segmentation1) || !.scalarString(segmentation2) ||
        !all(c(segmentation1, segmentation2) %in%
             names(object@segmentations))) {
        stop("segmentation selectors must be existing character names",
             call. = FALSE)
    }
    if (identical(segmentation1, segmentation2)) {
        stop("Compared segmentation names must be distinct", call. = FALSE)
    }
    methods::validObject(object)
    .validateUniverse(S4Vectors::SimpleList(molecules = object@molecules,
        features = object@features,
        fingerprint_value = object@universe_fingerprint,
        fingerprint_version = object@fingerprint_version,
        fingerprint_algorithm = object@fingerprint_algorithm))
    nm <- c(segmentation1, segmentation2)
    p <- object@segmentations[[segmentation1]]
    r <- object@segmentations[[segmentation2]]
    .checkPartitionIntegrity(p)
    .checkPartitionIntegrity(r)
    .compareAlignedPartitions(object, nm, p, r)
}
