# Assisted-by: OpenAI Codex
# Error-free floating expansions compare sums of the supplied binary doubles.
# A rounded scalar sum can conceal a real improvement smaller than its ULP.
.geometryExactAdd <- function(x, y) {
    out <- numeric()
    for (a in x) {
        total <- y + a
        recovered <- total - y
        error <- (y - (total - recovered)) + (a - recovered)
        if (error != 0) out <- c(out, error)
        y <- total
    }
    if (y != 0 || !length(out)) out <- c(out, y)
    out
}

.geometryExactSum <- function(x) {
    out <- 0
    for (value in x) out <- .geometryExactAdd(out, value)
    out
}

.geometryExactCompare <- function(x, y) {
    difference <- x
    for (value in y) difference <- .geometryExactAdd(difference, -value)
    sign(difference[length(difference)])
}

.geometryOptimum <- function(q, s, w) {
    if (!length(w)) return(integer())
    # Compact positional nodes, never scientific identifier conversion.
    q <- match(q, unique(q))
    s <- match(s, unique(s))
    n1 <- max(q)
    n2 <- max(s)
    selected <- rep(FALSE, length(w))
    repeat {
        distances <- vector("list", n1 + n2)
        previous <- rep(NA_integer_, n1 + n2)
        roots <- setdiff(seq_len(n1), q[selected])
        for (i in roots) distances[[i]] <- 0
        # Residual costs reverse selected edges; Bellman-Ford handles negatives.
        from <- ifelse(selected, n1 + s, q)
        to <- ifelse(selected, q, n1 + s)
        cost <- ifelse(selected, w, -w)
        for (iteration in seq_len(n1 + n2 - 1L)) {
            changed <- FALSE
            for (e in seq_along(w)) {
                if (is.null(distances[[from[e]]])) next
                candidate <- .geometryExactAdd(distances[[from[e]]], cost[e])
                old <- distances[[to[e]]]
                if (is.null(old) || .geometryExactCompare(candidate, old) < 0) {
                    distances[[to[e]]] <- candidate
                    previous[to[e]] <- e
                    changed <- TRUE
                }
            }
            if (!changed) break
        }
        terminals <- n1 + setdiff(seq_len(n2), s[selected])
        terminal <- NA_integer_
        best <- 0
        for (node in terminals) {
            value <- distances[[node]]
            if (!is.null(value) && .geometryExactCompare(value, best) < 0) {
                terminal <- node
                best <- value
            }
        }
        if (is.na(terminal)) break
        path <- integer()
        node <- terminal
        while (!is.na(previous[node])) {
            e <- previous[node]
            path <- c(path, e)
            node <- from[e]
            if (length(path) > n1 + n2) stop("Matching residual path cycle")
        }
        selected[path] <- !selected[path]
    }
    which(selected)
}

.geometryMatchingComponent <- function(q, s, w) {
    optimum <- .geometryOptimum(q, s, w)
    objective <- .geometryExactSum(w[optimum])
    chosen <- integer()
    remaining <- seq_along(w)
    for (e in seq_along(w)) {
        current <- .geometryExactSum(w[chosen])
        if (.geometryExactCompare(current, objective) == 0) {
            break
        }
        if (!e %in% remaining) next
        available <- remaining[remaining > e & q[remaining] != q[e] &
                                   s[remaining] != s[e]]
        residual <- .geometryOptimum(q[available], s[available], w[available])
        candidate <- c(chosen, e, available[residual])
        total <- .geometryExactSum(w[candidate])
        # Determinism requires EXACT equality of the unmodified objective.
        if (.geometryExactCompare(total, objective) == 0) {
            chosen <- c(chosen, e)
            remaining <- available
        } else {
            remaining <- remaining[remaining != e]
        }
    }
    if (.geometryExactCompare(.geometryExactSum(w[chosen]), objective) != 0) {
        stop("Exact optimal matching tie resolution failed")
    }
    chosen
}

.geometryMatching <- function(q, s, w) {
    if (!length(w)) return(integer())
    if (any(!is.finite(w) | w <= 0 | w > 1)) {
        stop("Matching requires finite positive IoU weights")
    }
    order <- order(q, s, method = "radix")
    q <- q[order]
    s <- s[order]
    w <- w[order]
    # Geometric components are computational partitions, not primary components.
    compact_q <- match(q, unique(q))
    compact_s <- match(s, unique(s))
    groups <- .graphComponents(compact_q, compact_s,
                              max(compact_q), max(compact_s))$side1[compact_q]
    result <- integer()
    for (group in unique(groups)) {
        index <- which(groups == group)
        selected <- .geometryMatchingComponent(q[index], s[index], w[index])
        result <- c(result, index[selected])
    }
    sort(order[result])
}
