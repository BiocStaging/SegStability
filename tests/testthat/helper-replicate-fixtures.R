wp10_set <- function(assignments = list(run_a = c(1L, 1L, 2L, 2L),
                                       run_b = c(1L, 1L, 2L, 2L)),
                     cell_ids = NULL, molecule_ids = NULL,
                     census_orders = NULL, molecule_orders = NULL,
                     run_metadata = NULL) {
    n <- length(assignments[[1L]])
    if (is.null(molecule_ids)) molecule_ids <- paste0("m", seq_len(n))
    if (!n) molecule_ids <- character()
    molecules <- S4Vectors::DataFrame(molecule_id = molecule_ids,
                                      feature_id = rep("g", n))
    runs <- lapply(seq_along(assignments), function(i) {
        a <- assignments[[i]]
        nm <- names(assignments)[i]
        ids <- if (is.null(cell_ids)) {
            paste0(nm, "_cell", seq_len(max(c(0L, a)) + 1L))
        } else cell_ids[[i]]
        census <- if (is.null(census_orders)) seq_along(ids) else census_orders[[i]]
        order <- if (is.null(molecule_orders)) seq_len(n) else molecule_orders[[i]]
        position <- a
        position[a <= 0L] <- 1L
        values <- ids[position]
        if (any(a <= 0L)) values[a <= 0L] <- NA
        status <- ifelse(a > 0L, "assigned",
                          ifelse(a == 0L, "unassigned", "ambiguous"))
        status <- as.character(status)
        metadata <- if (is.null(run_metadata)) S4Vectors::SimpleList() else {
            run_metadata[[i]]
        }
        SegmentationPartition(molecules, nm,
            S4Vectors::DataFrame(cell_id = ids[census]),
            molecule_ids[order], status[order], values[order],
            run_metadata = metadata)
    })
    names(runs) <- names(assignments)
    SegmentationSet(molecules, runs, "specimen")
}

# Independent tiny oracle: enumerate molecule pairs, not contingency margins.
wp10_pair_oracle <- function(left, right) {
    joint <- which(left > 0L & right > 0L)
    if (length(joint) < 2L) return(NA_real_)
    pairs <- utils::combn(joint, 2L)
    same_left <- left[pairs[1L, ]] == left[pairs[2L, ]]
    same_right <- right[pairs[1L, ]] == right[pairs[2L, ]]
    a <- sum(same_left & same_right)
    b <- sum(same_left & !same_right)
    c <- sum(!same_left & same_right)
    d <- sum(!same_left & !same_right)
    denominator <- (a + b) * (b + d) + (a + c) * (c + d)
    if (denominator == 0) NA_real_ else 2 * (a * d - b * c) / denominator
}

wp10_ll <- function(x) methods::as(as.character(x), "LLint")
