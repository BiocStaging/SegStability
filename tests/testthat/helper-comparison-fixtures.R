wp04_set <- function(matrix = base::matrix(c(3L, 0L, 0L, 2L), 2L),
                     extra = NULL, ids1 = NULL, ids2 = NULL) {
    n1 <- nrow(matrix)
    n2 <- ncol(matrix)
    if (is.null(ids1)) ids1 <- sprintf("x%d", seq_len(n1))
    if (is.null(ids2)) ids2 <- sprintf("a%d", seq_len(n2))
    # Reverse observation order deliberately; canonical edges must be sorted.
    pairs <- which(matrix > 0, arr.ind = TRUE)
    copies <- rep(seq_len(nrow(pairs)), matrix[pairs])
    i <- rev(pairs[copies, 1L])
    j <- rev(pairs[copies, 2L])
    s1 <- rep("assigned", length(i))
    s2 <- s1
    if (!is.null(extra)) {
        i <- c(i, extra$i)
        j <- c(j, extra$j)
        s1 <- c(s1, extra$s1)
        s2 <- c(s2, extra$s2)
    }
    m <- S4Vectors::DataFrame(molecule_id = sprintf("m%d", seq_along(i)),
                             feature_id = rep("g", length(i)))
    # Indexed extraction with NA is unsupported by LLint; fill only assigned.
    cells <- function(ids, index, status) {
        key <- rep(NA_character_, length(index))
        key[status == "assigned"] <- as.character(ids[index[
            status == "assigned"]])
        switch(class(ids)[1L], integer = as.integer(key), character = key,
               LLint = S4Vectors::as.LLint(key),
               integer64 = bit64::as.integer64(key))
    }
    p <- SegmentationPartition(m, "left",
        S4Vectors::DataFrame(cell_id = ids1), m$molecule_id, s1,
        cells(ids1, i, s1))
    q <- SegmentationPartition(m, "right",
        S4Vectors::DataFrame(cell_id = ids2), m$molecule_id, s2,
        cells(ids2, j, s2))
    SegmentationSet(m, list(left = p, right = q), "specimen")
}

wp04_compare <- function(...) compareSegmentations(wp04_set(...),
                                                  "left", "right")

wp04_class_input <- function() {
    h <- S4Vectors::Hits(1L, 1L, 2L, 2L)
    S4Vectors::mcols(h) <- S4Vectors::DataFrame(weight = 2,
        fraction_1 = 1, fraction_2 = 1, dominant_1 = TRUE,
        dominant_2 = TRUE, unique_dominant_1 = TRUE,
        unique_dominant_2 = TRUE, reciprocal_dominant = TRUE)
    cells <- S4Vectors::DataFrame(cell_id = c("used", "empty"),
        assigned_mass = c(2, 0), matched_mass = c(2, 0),
        unmatched_mass = c(0, 0), partner_count = c(1L, 0L),
        dominance_count = c(1L, 0L), component_id = c(1L, 2L),
        matched_fraction = c(1, NA_real_))
    other <- cells
    other$component_id <- c(1L, 3L)
    status <- rep(c("assigned", "unassigned", "ambiguous"), each = 3L)
    list(segmentation_names = c("left", "right"),
        cell_summary_1 = cells, cell_summary_2 = other, edges = h,
        components = S4Vectors::DataFrame(component_id = 1:3,
            n_segmentation_1 = c(1L, 1L, 0L),
            n_segmentation_2 = c(1L, 0L, 1L),
            component_class = c("one_to_one", "unmatched", "unmatched")),
        assignment_summary = S4Vectors::DataFrame(status_1 = status,
            status_2 = rep(c("assigned", "unassigned", "ambiguous"), 3L),
            n_molecules = c(2, rep(0, 8))), universe_size = 2,
        universe_fingerprint = strrep("a", 64L),
        assignment_fingerprints = c(strrep("b", 64L), strrep("c", 64L)),
        fingerprint_version = "1", fingerprint_algorithm = "sha256",
        schema_version = "1", coordinate_info = S4Vectors::SimpleList(),
        provenance = S4Vectors::SimpleList(), metadata = S4Vectors::SimpleList())
}
