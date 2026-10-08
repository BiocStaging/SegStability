test_that("empty coordinate declarations retain unknown semantics", {
    expect_silent(SegStability:::.adapterSpecimen(S4Vectors::SimpleList(), "roi"))
    expect_error(SegStability:::.adapterSpecimen(
        S4Vectors::SimpleList(specimen_id = "other"), "roi"), "conflict")
})

test_that("table ingestion preserves exact ID payloads and rejects lossy conversion", {
    for (kind in c("character", "integer", "integer64")) {
        if (kind == "integer64") skip_if_not_installed("bit64")
        ids <- wp11_ids(kind)
        tab <- data.frame(molecule_id = ids, feature_id = rep("g", 3L))
        expect_identical(SegStability:::.ingestionTable(tab, "table")$molecule_id, ids)
    }
    tab <- data.frame(a = 1L, b = 2L)
    names(tab) <- c("id", "id")
    expect_error(SegStability:::.ingestionTable(tab, "table"), "duplicate")
})

test_that("geometry and source adapters never invoke analytical entrypoints", {
    source <- lapply(c("fromMoleculeExperiment", "geometryFromSpatialFeatureExperiment"),
        function(name) paste(deparse(get(name, asNamespace("SegStability"))), collapse = "\n"))
    forbidden <- c("compareSegmentations", "stateStability", "supportSensitivity",
        "anchorTopology", "compareGeometry", "replicateAgreement")
    for (name in forbidden) expect_false(any(grepl(paste0(name, "\\("), source)))
})

test_that("base table ingestion retains LLint IDs including typed missing cell entries", {
    ids <- wp11_ids("LLint")
    molecules <- structure(list(molecule_id = ids, feature_id = rep("g", 3L)),
        class = "data.frame", row.names = c(NA_integer_, -3L))
    cells <- structure(list(cell_id = ids), class = "data.frame",
        row.names = c(NA_integer_, -3L))
    p <- SegmentationPartition(molecules, "native", cells, ids,
        c("assigned", "unassigned", "ambiguous"),
        wp11_subset(ids, c(1L, NA_integer_, NA_integer_)))
    expect_identical(cellCensus(p)$cell_id, ids)
    ss <- SegmentationSet(molecules, list(native = p), "roi")
    expect_identical(moleculeData(ss)$molecule_id, ids)
})

test_that("retained specimen declarations are checked without inference", {
    check <- SegStability:::.adapterSourceSpecimen
    expect_silent(check(S4Vectors::DataFrame(specimen_id = c("roi", NA_character_)), "roi"))
    expect_silent(check(S4Vectors::DataFrame(note = "opaque"), "roi"))
    expect_error(check(S4Vectors::DataFrame(specimen_id = "other"), "roi"), "specimen")
    expect_error(check(S4Vectors::DataFrame(specimen_id = 1), "roi"), "specimen")
})

test_that("explicit universe linking preserves every scientific ID family", {
    for (kind in c("character", "integer", "LLint", "integer64")) {
        if (kind == "integer64") skip_if_not_installed("bit64")
        source <- wp11_tables(kind)$molecules
        u <- source[c(3L, 1L, 2L), ]
        expect_identical(SegStability:::.meUniverse(source, u), u)
        extra <- source
        extra$molecule_id <- source$molecule_id[c(1L, 1L, 3L)]
        expect_error(SegStability:::.meUniverse(extra, u), "Duplicate")
        u$feature_id[1L] <- "conflicting feature"
        expect_error(SegStability:::.meUniverse(source, u), "feature")
    }
})
