wp08_info <- function() S4Vectors::SimpleList(units = "micrometres",
    origin = c(0, 0), axis_directions = c("right", "up"),
    transforms = S4Vectors::SimpleList(), registration = "common local XY",
    specimen_id = "specimen")

wp08_window <- function() c(xmin = 0, xmax = 10, ymin = 0, ymax = 10)

wp08_box <- function(xmin = 1, xmax = 2, ymin = 1, ymax = 2) {
    sf::st_polygon(list(rbind(c(xmin, ymin), c(xmax, ymin), c(xmax, ymax),
                             c(xmin, ymax), c(xmin, ymin))))
}

wp08_geometry <- function(cmp) {
    out <- lapply(segmentationNames(cmp), function(side) {
        ids <- cellSummary(cmp, side)$cell_id
        shapes <- lapply(seq_along(ids), function(i) wp08_box())
        g <- sf::st_sf(geometry = sf::st_sfc(shapes))
        # Preserve typed identifiers without sf constructor coercion.
        g$cell_id <- ids
        g
    })
    names(out) <- segmentationNames(cmp)
    out
}

wp08_run <- function(matrix, ...) {
    cmp <- wp04_compare(matrix, ...)
    anchorTopology(cmp, wp08_geometry(cmp), wp08_window(), wp08_info())
}

wp08_unchanged <- function(cmp, g, expr) {
    before <- serialize(list(cmp, g), NULL)
    force(expr)
    expect_identical(serialize(list(cmp, g), NULL), before)
}
