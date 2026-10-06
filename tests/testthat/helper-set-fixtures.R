wp03_set <- function(...) {
    run <- S4Vectors::SimpleList(seed = 7L, nested = S4Vectors::SimpleList(
        list(parameters = list(a = 1L, b = c("x", "y")))),
        table = S4Vectors::DataFrame(label = c("a", NA_character_)),
        unknown = NA_real_, bytes = as.raw(c(0, 255)))
    args <- list(molecules = wp02_molecules(),
        segmentations = list(native = wp02_partition(run_metadata = run),
            other = wp02_partition(name = "other",
                cells = S4Vectors::DataFrame(cell_id = c("c1", "extra")))),
        specimen_id = "specimen")
    changes <- list(...)
    args[names(changes)] <- changes
    do.call(SegmentationSet, args)
}
