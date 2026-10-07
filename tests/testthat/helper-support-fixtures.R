wp07_state_fields <- function() {
    c("resolved_mass", "same_state_mass", "different_state_mass",
      "unresolved_label_mass", "conditional_stability",
      "same_state_fraction_of_matched_mass", "total_mass_stability",
      "endpoint_unavailable_cells", "total_estimable_cells",
      "total_nonestimable_cells", "conditional_estimable_cells",
      "conditional_nonestimable_cells", "mean_cell_total_mass_stability",
      "mean_cell_conditional_stability")
}

wp07_conservation <- function(out, i, cmp) {
    r <- out[i, ]
    expect_identical(r$retained_edges + r$removed_edges, r$original_edges)
    expect_identical(r$retained_mass + r$excluded_matched_mass, r$original_matched_mass)
    components <- r$components[[1]]
    expect_identical(sum(components$selected_retained_mass), r$retained_mass)
    expect_identical(sum(components$selected_removed_mass), r$excluded_matched_mass)
    for (side in 1:2) {
        c <- r$cells[[1]][[side]]
        original <- cellSummary(cmp, segmentationNames(cmp)[side])
        expect_identical(c[, colnames(original)], original)
        expect_identical(sum(c$selected_retained_mass), r$retained_mass)
        expect_identical(sum(c$selected_removed_mass), r$excluded_matched_mass)
        expect_identical(c$selected_retained_mass + c$selected_removed_mass, c$matched_mass)
        expect_true(all(c$selected_retained_mass <= c$assigned_mass))
    }
    if (r$states_supplied) {
        expect_identical(r$same_state_mass + r$different_state_mass +
                             r$unresolved_label_mass, r$retained_mass)
        expect_identical(sum(components$selected_resolved_mass), r$resolved_mass)
        for (side in 1:2) {
            expect_identical(sum(r$cells[[1]][[side]]$selected_same_state_mass),
                             r$same_state_mass)
        }
    }
}
