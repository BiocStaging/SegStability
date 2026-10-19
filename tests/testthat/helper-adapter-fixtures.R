# Assisted-by: OpenAI Codex
wp11_ids <- function(kind = "character") {
    text <- c("9007199254740992", "9007199254740993", "9223372036854775807")
    switch(kind, character = c("001", " 1 ", "e\u0301"), integer = 1:3,
        LLint = S4Vectors::as.LLint(text), integer64 = bit64::as.integer64(text))
}

wp11_subset <- function(ids, index) {
    if (methods::is(ids, "LLint")) {
        return(S4Vectors::as.LLint(as.character(ids)[index]))
    }
    ids[index]
}

wp11_tables <- function(kind = "character") {
    ids <- wp11_ids(kind)
    list(molecules = S4Vectors::DataFrame(molecule_id = ids,
            feature_id = c(" G ", " G ", "g"),
            x = c(1, 2, 3), y = c(4, 5, 6)),
        cells = S4Vectors::DataFrame(cell_id = ids, note = c("a", "zero", "zero")),
        assignments = S4Vectors::DataFrame(molecule_id = ids[c(3L, 1L, 2L)],
            status = c("ambiguous", "assigned", "unassigned"),
            cell_id = wp11_subset(ids, c(NA_integer_, 1L, NA_integer_))))
}

wp11_me <- function(tables = wp11_tables()) {
    skip_if_not_installed("MoleculeExperiment")
    m <- tables$molecules
    leaves <- lapply(unique(m$feature_id), function(feature) {
        out <- m[m$feature_id == feature, , drop = FALSE]
        out$feature_id <- NULL
        names(out)[names(out) == "x"] <- "x_location"
        names(out)[names(out) == "y"] <- "y_location"
        out
    })
    names(leaves) <- unique(m$feature_id)
    MoleculeExperiment::MoleculeExperiment(
        molecules = list(detected = list(roi = leaves)))
}

wp11_import <- function(object = wp11_me(), tables = wp11_tables(), ...) {
    args <- list(object = object, assay = "detected", sample_id = "roi",
        name = "native", specimen_id = "specimen", cells = tables$cells,
        assignments = tables$assignments)
    changes <- list(...)
    args[names(changes)] <- changes
    do.call(fromMoleculeExperiment, args)
}

wp11_manual <- function(tables, ambiguity = NULL,
                        run_metadata = S4Vectors::SimpleList(),
                        provenance = S4Vectors::SimpleList(),
                        metadata = S4Vectors::SimpleList(),
                        coordinate_info = S4Vectors::SimpleList()) {
    a <- tables$assignments
    p <- SegmentationPartition(tables$molecules, "native", tables$cells,
        a$molecule_id, a$status, a$cell_id, ambiguity, run_metadata,
        provenance, metadata)
    SegmentationSet(tables$molecules, list(native = p), "specimen",
        coordinate_info, provenance, metadata)
}

wp11_sfe <- function(kind = "character", invalid = FALSE) {
    skip_if_not_installed("SpatialFeatureExperiment")
    skip_if_not_installed("sf")
    ring <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
    if (invalid) ring <- rbind(c(0, 0), c(1, 1), c(0, 1), c(1, 0), c(0, 0))
    g <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(ring)),
        sf::st_polygon(list(ring + 2))))
    rownames(g) <- c("source-A", "source-B")
    sfe <- SpatialFeatureExperiment::SpatialFeatureExperiment(
        assays = list(counts = matrix(0L, 1L, 2L,
            dimnames = list("g", rownames(g)))),
        colData = S4Vectors::DataFrame(sample_id = rep("roi", 2L)),
        spatialCoords = matrix(c(0.5, 2.5, 0.5, 2.5), 2L,
            dimnames = list(NULL, c("x", "y"))),
        sample_id = "roi", unit = "micron")
    SpatialFeatureExperiment::colGeometry(sfe, "cellSeg", sample_id = "roi",
        translate = FALSE) <- g
    sfe
}

wp11_geometry <- function(object = wp11_sfe(), kind = "character", ...) {
    ids <- wp11_ids(kind)
    args <- list(object = object, type = "cellSeg", sample_id = "roi",
        geometry_kind = "col", cell_map = S4Vectors::DataFrame(
            geometry_id = c("source-B", "source-A"), cell_id = ids[c(2L, 1L)]),
        cells = S4Vectors::DataFrame(cell_id = ids),
        coordinate_info = S4Vectors::SimpleList(units = "micron",
            specimen_id = "specimen"))
    changes <- list(...)
    args[names(changes)] <- changes
    do.call(geometryFromSpatialFeatureExperiment, args)
}
