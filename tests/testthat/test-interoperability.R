# Assisted-by: OpenAI Codex
test_that("the approved adapter functions exist with exact interfaces", {
    expect_true(is.function(fromMoleculeExperiment))
    expect_identical(names(formals(fromMoleculeExperiment)), c("object", "assay",
        "sample_id", "name", "specimen_id", "cells", "assignments", "molecule_id_col",
        "coordinate_cols", "universe", "ambiguity", "coordinate_info", "run_metadata",
        "provenance", "metadata"))
    expect_true(is.function(geometryFromSpatialFeatureExperiment))
    expect_identical(names(formals(geometryFromSpatialFeatureExperiment)), c("object",
        "type", "sample_id", "geometry_kind", "cell_map", "cells", "coordinate_info",
        "reference"))
})

test_that("base data frames preserve scientific constructor semantics", {
    t <- wp11_tables()
    manual <- wp11_manual(t)
    t$molecules <- as.data.frame(t$molecules)
    t$cells <- as.data.frame(t$cells)
    expect_identical(wp11_manual(t), manual)
})

test_that("ME ingestion exactly equals manual construction in every ID family", {
    for (kind in c("character", "integer", "LLint", "integer64")) {
        if (kind == "integer64") skip_if_not_installed("bit64")
        t <- wp11_tables(kind)
        me <- wp11_me(t)
        before <- serialize(me, NULL)
        args <- wp11_import(me, t)
        expect_identical(names(args), c("molecules", "segmentations", "specimen_id",
            "coordinate_info", "provenance", "metadata"))
        set <- do.call(SegmentationSet, as.list(args))
        expect_identical(set, wp11_manual(t))
        expect_identical(moleculeData(set)$molecule_id, t$molecules$molecule_id)
        expect_identical(cellCensus(args$segmentations[[1L]]), t$cells)
        expect_identical(assignmentStatus(args$segmentations[[1L]]),
            c("assigned", "unassigned", "ambiguous"))
        expect_identical(serialize(me, NULL), before)
        path <- tempfile(fileext = ".rds")
        saveRDS(set, path)
        expect_identical(readRDS(path), set)
        unlink(path)
    }
})

test_that("ME ingestion preserves ambiguity and supplied metadata", {
    t <- wp11_tables()
    ambiguity <- S4Vectors::SimpleList(evidence = S4Vectors::DataFrame(
        molecule_id = t$molecules$molecule_id[3L], candidate_count = 2L,
        candidate_ids_known = TRUE), candidates = S4Vectors::DataFrame(
        molecule_id = rep(t$molecules$molecule_id[3L], 2L), cell_id = t$cells$cell_id[2:3]))
    run <- S4Vectors::SimpleList(seed = 7L)
    prov <- S4Vectors::SimpleList(software = "explicit", software_version = "1")
    meta <- S4Vectors::SimpleList(note = "opaque")
    args <- wp11_import(ambiguity = ambiguity, run_metadata = run,
        provenance = prov, metadata = meta)
    expect_identical(do.call(SegmentationSet, as.list(args)),
        wp11_manual(t, ambiguity, run, prov, meta))
    for (count in list(2L, NA_integer_)) {
        ambiguity$evidence$candidate_count <- count
        ambiguity$evidence$candidate_ids_known <- FALSE
        ambiguity$candidates <- ambiguity$candidates[FALSE, ]
        args <- wp11_import(ambiguity = ambiguity)
        expect_identical(do.call(SegmentationSet, as.list(args)), wp11_manual(t, ambiguity))
    }
})

test_that("ME fails when original identity or semantic inputs are absent", {
    me <- wp11_me()
    expect_error(wp11_import(me, molecule_id_col = "lost_id"), "molecule.*ID|molecule_id")
    expect_error(wp11_import(me, assay = "absent"), "assay|Assay")
    expect_error(wp11_import(me, sample_id = "absent"), "sample")
    expect_error(wp11_import(me, assignments = NULL), "assignments")
    expect_error(wp11_import(me, cells = NULL), "cells|census")
    expect_error(wp11_import(me, specimen_id = ""), "specimen")
    t <- wp11_tables()
    t$assignments$status[1L] <- "low_confidence"
    expect_error(wp11_import(me, t), "status|Status")
    t <- wp11_tables()
    t$assignments$cell_id[1L] <- t$cells$cell_id[1L]
    expect_error(wp11_import(me, t), "unresolved|ambiguous|missing")
})

test_that("ME aligns an explicitly supplied universe without guessing", {
    t <- wp11_tables()
    u <- t$molecules[c(3L, 1L, 2L), ]
    out <- wp11_import(universe = u)
    t$molecules <- u
    expect_identical(do.call(SegmentationSet, as.list(out)), wp11_manual(t))
    u$feature_id[1L] <- "conflict"
    expect_error(wp11_import(universe = u), "feature")
    u <- t$molecules
    u$x[1L] <- 99
    expect_error(wp11_import(universe = u), "coordinate")
    expect_error(wp11_import(universe = u[-1L, ]), "universe|missing|extra")
})

test_that("ME rejects unsafe ID classes and malformed coordinate mapping", {
    for (bad in list(c(1, 2, 3), factor(c("1", "2", "3")), c(TRUE, FALSE, TRUE))) {
        t <- wp11_tables()
        t$molecules$molecule_id <- bad
        expect_error(wp11_import(wp11_me(t), t), "ID|representation")
    }
    t <- wp11_tables()
    t$molecules$x[1L] <- Inf
    expect_error(wp11_import(wp11_me(t), t), "finite|coordinate")
    expect_error(wp11_import(coordinate_cols = c(x = "x_location")), "coordinate")
    expect_error(wp11_import(coordinate_info = S4Vectors::SimpleList(
        specimen_id = "other")), "specimen|conflict")
    out <- wp11_import(coordinate_cols = NULL)
    expect_false(any(c("x", "y", "z") %in% names(out$molecules)))
})

test_that("SFE sidecars preserve keyed geometry and exact cell IDs", {
    for (kind in c("character", "integer", "LLint", "integer64")) {
        if (kind == "integer64") skip_if_not_installed("bit64")
        sfe <- wp11_sfe(kind)
        before <- serialize(sfe, NULL)
        original <- SpatialFeatureExperiment::colGeometry(sfe, "cellSeg", sample_id = "roi")
        out <- wp11_geometry(sfe, kind)
        expect_identical(names(out), c("geometry", "coordinate_info", "validation"))
        expect_identical(out$geometry$cell_id, wp11_ids(kind)[1:2])
        expect_identical(sf::st_geometry(out$geometry), sf::st_geometry(original))
        expect_identical(sf::st_crs(out$geometry), sf::st_crs(original))
        expect_true(any(out$validation$status == "UNKNOWN"))
        expect_false(any(out$validation$status == "FAIL" & out$validation$severity == "error"))
        expect_identical(serialize(sfe, NULL), before)
    }
})

test_that("SFE reports invalid geometry and rejects identity/mapping defects", {
    out <- wp11_geometry(wp11_sfe(invalid = TRUE))
    expect_true(any(out$validation$status == "FAIL"))
    expect_identical(nrow(out$geometry), 2L)
    expect_error(wp11_geometry(geometry_kind = "row"), "geometry_kind")
    expect_error(wp11_geometry(type = "absent"), "geometry|type|name|Name")
    expect_error(wp11_geometry(cell_map = S4Vectors::DataFrame(
        geometry_id = "source-A", cell_id = "001")), "map|key|missing")
    expect_error(wp11_geometry(cell_map = S4Vectors::DataFrame(
        geometry_id = c("source-A", "source-A"), cell_id = c("001", " 1 "))), "Duplicate|duplicate")
    expect_error(wp11_geometry(cell_map = S4Vectors::DataFrame(
        geometry_id = c("source-A", "source-B"), cell_id = c(1, 2))), "ID|representation")
})

test_that("unsupported objects and absent packages fail informatively", {
    expect_error(fromMoleculeExperiment(NULL), "MoleculeExperiment")
    expect_error(geometryFromSpatialFeatureExperiment(NULL), "SpatialFeatureExperiment")
})

test_that("ME preserves ancillary tables, rejects conflicting specimen facts and no ID repair", {
    t <- wp11_tables()
    t$molecules$specimen_id <- rep("other", 3L)
    expect_error(wp11_import(wp11_me(t), t), "specimen")
    t <- wp11_tables()
    t$molecules$note <- c(" Keep ", "KEEP", "keep")
    out <- wp11_import(wp11_me(t), t)
    expect_identical(out$molecules$note, t$molecules$note)
    expect_identical(do.call(SegmentationSet, as.list(out)), wp11_manual(t))
    t$assignments$molecule_id[1L] <- "\u00e9"
    expect_error(wp11_import(wp11_me(t), t), "mismatch|missing|extra")
})

test_that("ME supports retained z and explicit typed empty tables", {
    t <- wp11_tables()
    t$molecules$z_location <- c(7, 8, 9)
    args <- wp11_import(wp11_me(t), t, coordinate_cols =
        c(x = "x_location", y = "y_location", z = "z_location"))
    expect_identical(args$molecules$z, t$molecules$z_location)
    for (kind in c("character", "integer", "LLint", "integer64")) {
        if (kind == "integer64") skip_if_not_installed("bit64")
        t <- wp11_tables(kind)
        t$molecules <- t$molecules[FALSE, ]
        t$assignments <- t$assignments[FALSE, ]
        leaf <- S4Vectors::DataFrame(molecule_id = wp11_ids(kind)[FALSE],
            x_location = numeric(), y_location = numeric())
        skip_if_not_installed("MoleculeExperiment")
        me <- MoleculeExperiment::MoleculeExperiment(
            molecules = list(detected = list(roi = list(g = leaf))))
        args <- wp11_import(me, t)
        expect_identical(do.call(SegmentationSet, as.list(args)), wp11_manual(t))
        expect_identical(nrow(cellCensus(args$segmentations[[1L]])), 3L)
    }
})

test_that("SFE annotation geometry uses explicit retained keys without fallback", {
    sfe <- wp11_sfe()
    g <- SpatialFeatureExperiment::colGeometry(sfe, "cellSeg", sample_id = "roi")
    SpatialFeatureExperiment::annotGeometry(sfe, "cells", sample_id = "roi",
        translate = FALSE) <- g
    out <- wp11_geometry(sfe, type = "cells", geometry_kind = "annot")
    expect_identical(sf::st_geometry(out$geometry), sf::st_geometry(g))
    expect_identical(out$geometry$cell_id, wp11_ids()[1:2])
    expect_error(wp11_geometry(sfe, type = "cells", geometry_kind = "col"),
        "geometry|name|Name|type")
})

test_that("SFE preserves CRS and refuses conflicting coordinate declarations", {
    sfe <- wp11_sfe()
    g <- SpatialFeatureExperiment::colGeometry(sfe, "cellSeg", sample_id = "roi")
    sf::st_crs(g) <- 3857
    SpatialFeatureExperiment::colGeometry(sfe, "projected", sample_id = "roi",
        translate = FALSE) <- g
    out <- wp11_geometry(sfe, type = "projected", coordinate_info =
        S4Vectors::SimpleList(units = "metre", specimen_id = "specimen", crs = "EPSG:3857"))
    expect_identical(sf::st_crs(out$geometry), sf::st_crs(g))
    expect_identical(sf::st_geometry(out$geometry), sf::st_geometry(g))
    expect_error(wp11_geometry(sfe, type = "projected", coordinate_info =
        S4Vectors::SimpleList(units = "metre", specimen_id = "specimen", crs = "EPSG:4326")),
        "CRS|crs|conflict")
    reference <- wp11_manual(wp11_tables())
    expect_error(wp11_geometry(sfe, reference = reference, coordinate_info =
        S4Vectors::SimpleList(units = "micron", specimen_id = "other")), "conflict")
})

test_that("ME rejects heterogeneous tables instead of coercing or silently dropping records", {
    skip_if_not_installed("MoleculeExperiment")
    leaf <- S4Vectors::DataFrame(molecule_id = "a", x_location = 1, y_location = 2)
    other <- leaf
    other$molecule_id <- 1L
    me <- MoleculeExperiment::MoleculeExperiment(
        molecules = list(detected = list(roi = list(g = leaf, h = other))))
    expect_error(wp11_import(me), "famil")
    other <- leaf
    other$molecule_id <- "b"
    other$extra <- "metadata"
    me <- MoleculeExperiment::MoleculeExperiment(
        molecules = list(detected = list(roi = list(g = leaf, h = other))))
    expect_error(wp11_import(me), "schema")
    other <- leaf
    me <- MoleculeExperiment::MoleculeExperiment(
        molecules = list(detected = list(roi = list(g = leaf, h = other))))
    expect_error(wp11_import(me), "Duplicate")
})

test_that("ME retains character IDs outside signed 64-bit range without normalization", {
    t <- wp11_tables()
    ids <- c("9223372036854775809", "9223372036854775810", "0009223372036854775809")
    t$molecules$molecule_id <- ids
    t$cells$cell_id <- ids
    t$assignments$molecule_id <- ids[c(3L, 1L, 2L)]
    t$assignments$cell_id <- ids[c(NA_integer_, 1L, NA_integer_)]
    out <- wp11_import(wp11_me(t), t)
    expect_identical(do.call(SegmentationSet, as.list(out)), wp11_manual(t))
    expect_identical(out$molecules$molecule_id, ids)
})

test_that("cell/spot containers are not molecule adapters and unnamed geometry is unsupported", {
    sfe <- wp11_sfe()
    spatial <- methods::as(sfe, "SpatialExperiment")
    expect_error(wp11_import(spatial), "Expected a MoleculeExperiment")
    expect_error(wp11_import(sfe), "Expected a MoleculeExperiment")
    expect_error(wp11_geometry(spatial), "Expected a SpatialFeatureExperiment")
    colnames(sfe) <- NULL
    expect_error(wp11_geometry(sfe), "explicit source column keys")
})
