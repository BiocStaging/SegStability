wp06_states <- function(cmp, left, right) {
    maps <- list(left, right)
    for (i in seq_along(maps)) {
        maps[[i]] <- S4Vectors::DataFrame(
            cell_id = cellSummary(cmp, segmentationNames(cmp)[i])$cell_id,
            state = maps[[i]])
    }
    names(maps) <- segmentationNames(cmp)
    maps
}

wp06_result <- function(matrix, left, right, ...) {
    cmp <- wp04_compare(matrix, ...)
    stateStability(cmp, wp06_states(cmp, left, right))
}

wp06_conservation <- function(result) {
    g <- result$global
    expect_identical(g$same_state_mass + g$different_state_mass +
                         g$unresolved_label_mass, g$matched_mass)
    for (field in c("matched_mass", "resolved_mass", "same_state_mass",
                    "different_state_mass", "unresolved_label_mass")) {
        expect_identical(sum(result$components[[field]]), g[[field]])
        for (side in 1:2) {
            expect_identical(sum(result$cells[[side]][[field]]), g[[field]])
        }
    }
    for (side in 1:2) {
        c <- result$cells[[side]]
        expect_true(all(c$same_state_mass <= c$resolved_mass))
        expect_true(all(c$resolved_mass <= c$matched_mass))
        expect_true(all(c$matched_mass <= c$assigned_mass))
    }
}
