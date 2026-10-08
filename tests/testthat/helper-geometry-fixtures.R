wp09_geometry <- function(cmp, left, right) {
    shapes <- list(left, right)
    out <- lapply(1:2, function(i) {
        g <- sf::st_sf(geometry = sf::st_sfc(shapes[[i]]))
        g$cell_id <- cellSummary(cmp, segmentationNames(cmp)[i])$cell_id
        g
    })
    names(out) <- segmentationNames(cmp)
    out
}

wp09_run <- function(cmp, g) compareGeometry(cmp, g, wp08_info())

# Independent enumeration of every matching, including unmatched endpoints.
wp09_oracle <- function(q, s, w) {
    candidates <- list(integer())
    for (e in seq_along(w)) {
        additions <- lapply(candidates, function(x) {
            if (q[e] %in% q[x] || s[e] %in% s[x]) return(NULL)
            c(x, e)
        })
        candidates <- c(candidates, Filter(Negate(is.null), additions))
    }
    scores <- vapply(candidates, function(x) sum(w[x]), numeric(1))
    best <- candidates[scores == max(scores)]
    # Inputs are canonical pairs; lexicographic incidence prefers earlier edges.
    for (e in seq_along(w)) {
        included <- vapply(best, function(x) e %in% x, logical(1))
        if (any(included)) best <- best[included]
    }
    best[[1L]]
}
