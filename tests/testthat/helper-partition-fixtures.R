wp02_molecules <- function() S4Vectors::DataFrame(
    molecule_id = c("m1", "m2", "m3"), feature_id = c("g1", "g2", "g1"))

wp02_partition <- function(...) {
    args <- list(molecules = wp02_molecules(), name = "native",
        cells = S4Vectors::DataFrame(cell_id = c("c1", "c0", "unused")),
        molecule_id = c("m1", "m2", "m3"),
        status = c("assigned", "unassigned", "ambiguous"),
        cell_id = c("c1", NA_character_, NA_character_))
    changes <- list(...)
    args[names(changes)] <- changes
    do.call(SegmentationPartition, args)
}

wp02_ambiguity <- function() S4Vectors::SimpleList(
    evidence = S4Vectors::DataFrame(molecule_id = "m3",
        candidate_count = 2L, candidate_ids_known = TRUE),
    candidates = S4Vectors::DataFrame(molecule_id = c("m3", "m3"),
        cell_id = c("c0", "c1")))
