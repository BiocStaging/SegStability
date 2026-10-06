.canonicalStatus <- function(status) {
    if (!identical(class(status), "character") || anyNA(status) ||
        any(!status %in% c("assigned", "unassigned", "ambiguous"))) {
        stop("status must contain assigned, unassigned or ambiguous names",
             call. = FALSE)
    }
    base::match(status, c("assigned", "unassigned", "ambiguous"))
}

.mapCells <- function(cell_id, census, status) {
    .validateIds(census, "cell census")
    assigned <- status == 1L
    if (is.null(cell_id)) {
        if (any(assigned)) {
            stop("assigned molecules require cell_id", call. = FALSE)
        }
        return(rep(NA_integer_, length(status)))
    }
    if (length(cell_id) != length(status)) {
        stop("cell_id and status lengths must agree", call. = FALSE)
    }
    missing <- is.na(cell_id)
    if (any(missing[assigned])) {
        stop("assigned molecules require nonmissing cell_id", call. = FALSE)
    }
    if (any(!missing[!assigned])) {
        stop("unassigned/ambiguous molecules cannot have forced cell_id",
             call. = FALSE)
    }
    if (.idFamily(cell_id) != .idFamily(census)) {
        stop("cell_id and census ID families must agree", call. = FALSE)
    }
    out <- rep(NA_integer_, length(status))
    out[assigned] <- .matchIds(cell_id[which(assigned)], census)
    if (anyNA(out[assigned])) {
        stop("Unknown assigned cell_id outside census", call. = FALSE)
    }
    out
}

.candidateCounts <- function(x) {
    if (!(is.integer(x) || (is.double(x) && !is.object(x))) ||
        any(!is.na(x) & (!is.finite(x) | x < 0 | x != floor(x) |
                         x > .Machine$integer.max)) || any(is.nan(x))) {
        stop("candidate_count must contain nonnegative whole counts or NA",
             call. = FALSE)
    }
    as.integer(x)
}

.prepareAmbiguity <- function(ambiguity, molecules, cells, status) {
    ambiguous <- which(status == 3L)
    evidence <- S4Vectors::DataFrame(molecule_index = ambiguous,
        candidate_count = rep(NA_integer_, length(ambiguous)),
        candidate_ids_known = rep(FALSE, length(ambiguous)))
    from <- to <- integer()
    if (!is.null(ambiguity)) {
        if (!methods::is(ambiguity, "SimpleList") ||
            !identical(sort(names(ambiguity)), c("candidates", "evidence"))) {
            stop("ambiguity must contain evidence and candidates",
                 call. = FALSE)
        }
        e <- ambiguity$evidence
        c <- ambiguity$candidates
        if (!methods::is(e, "DataFrame") ||
            !all(c("molecule_id", "candidate_count", "candidate_ids_known")
                 %in% colnames(e)) || !methods::is(c, "DataFrame") ||
            !all(c("molecule_id", "cell_id") %in% colnames(c))) {
            stop("Invalid external ambiguity tables", call. = FALSE)
        }
        index <- .matchIds(e$molecule_id, molecules)
        if (anyDuplicated(index) || anyNA(index) ||
            any(!index %in% ambiguous)) {
            stop("Evidence requires unique ambiguous molecule references",
                 call. = FALSE)
        }
        counts <- .candidateCounts(e$candidate_count)
        if (!is.logical(e$candidate_ids_known) ||
            anyNA(e$candidate_ids_known)) {
            stop("candidate_ids_known must be nonmissing logical",
                 call. = FALSE)
        }
        rows <- base::match(index, ambiguous)
        evidence$candidate_count[rows] <- counts
        evidence$candidate_ids_known[rows] <- e$candidate_ids_known
        from <- .matchIds(c$molecule_id, molecules)
        to <- .matchIds(c$cell_id, cells)
        if (anyNA(from) || anyNA(to)) {
            stop("Unknown ambiguity candidate reference", call. = FALSE)
        }
        if (any(!from %in% ambiguous)) {
            stop("Candidate relationships require ambiguous molecules",
                 call. = FALSE)
        }
        if (anyDuplicated(paste(from, to, sep = ":"))) {
            stop("Duplicate ambiguity candidate relationships", call. = FALSE)
        }
        order <- base::order(from, to, method = "radix")
        from <- from[order]
        to <- to[order]
    }
    hits <- S4Vectors::Hits(from = from, to = to, nLnode = length(molecules),
                            nRnode = length(cells))
    observed <- tabulate(from, nbins = length(molecules))[ambiguous]
    known <- evidence$candidate_ids_known
    if (any(observed[!known] != 0L)) {
        stop("Candidates supplied while identities are declared unknown",
             call. = FALSE)
    }
    specified <- known & !is.na(evidence$candidate_count)
    if (any(evidence$candidate_count[specified] != observed[specified])) {
        stop("candidate_count disagrees with known candidates", call. = FALSE)
    }
    evidence$candidate_count[known] <- observed[known]
    S4Vectors::SimpleList(evidence = evidence, candidates = hits)
}
