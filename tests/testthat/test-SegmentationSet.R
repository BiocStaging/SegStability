test_that("set owns the frozen universe and preserves independent partitions", {
    x <- wp03_set()
    expect_s4_class(x, "SegmentationSet")
    expect_true(methods::validObject(x))
    expect_identical(sort(methods::slotNames(x)), sort(c("molecules", "features",
        "segmentations", "specimen_id", "coordinate_info", "provenance",
        "universe_fingerprint", "fingerprint_version", "fingerprint_algorithm",
        "schema_version", "metadata")))
    expect_s4_class(x@segmentations, "SimpleList")
    expect_identical(segmentationNames(x), c("native", "other"))
    expect_identical(moleculeData(x), wp02_molecules())
    expect_identical(featureData(x)$feature_id, c("g1", "g2"))
    expect_identical(cellCensus(x@segmentations$native)$cell_id,
                     c("c1", "c0", "unused"))
    expect_identical(cellCensus(x@segmentations$other)$cell_id,
                     c("c1", "extra"))
    expect_identical(assignmentStatus(x@segmentations$other),
                     c("assigned", "unassigned", "ambiguous"))
    expect_identical(x@fingerprint_version, "1")
    expect_identical(x@fingerprint_algorithm, "sha256")
    expect_identical(x@schema_version, "1")
    expect_false(any(c("molecules", "molecule_id", "features") %in%
        methods::slotNames(x@segmentations$native)))
    one <- wp03_set(segmentations = list(native = wp02_partition()))
    expect_identical(segmentationNames(one), "native")
    reversed <- wp03_set(segmentations = rev(x@segmentations))
    expect_identical(segmentationNames(reversed), c("other", "native"))
    expect_identical(reversed@universe_fingerprint, x@universe_fingerprint)
})

test_that("set constructor enforces exact scalar identity and collection names", {
    for (id in list(character(), c("a", "b"), NA_character_, "", 1L)) {
        expect_error(wp03_set(specimen_id = id), "specimen_id")
    }
    for (p in list(list(), S4Vectors::SimpleList(), list(wp02_partition()),
        setNames(list(wp02_partition()), ""),
        setNames(list(wp02_partition()), NA_character_),
        list(native = wp02_partition(), native = wp02_partition()),
        list(alias = wp02_partition()), list(native = 1L))) {
        expect_error(wp03_set(segmentations = p))
    }
    expect_error(wp03_set(segmentations = wp02_partition()), "list")
    for (field in c("coordinate_info", "metadata", "provenance")) {
        args <- setNames(list(list()), field)
        expect_error(do.call(wp03_set, args), "SimpleList")
    }
})

test_that("set rejects changed universe identity, feature map and assignment binding", {
    m <- wp02_molecules()
    for (changed in list(
        S4Vectors::DataFrame(molecule_id = c("m1", "m2", "other"),
            feature_id = m$feature_id), m[3:1, ],
        S4Vectors::DataFrame(molecule_id = 1:3, feature_id = m$feature_id),
        S4Vectors::DataFrame(molecule_id = m$molecule_id,
            feature_id = c("g2", "g1", "g2")), m[1:2, ])) {
        expect_error(wp03_set(molecules = changed), "binding")
    }
    bad <- wp02_partition()
    bad@cell_index[1L] <- 2L
    expect_true(methods::validObject(bad))
    expect_error(wp03_set(segmentations = list(native = bad)),
                 "assignment fingerprint")
    bad <- wp02_partition()
    bad@universe_fingerprint <- paste(rep("0", 64L), collapse = "")
    expect_error(wp03_set(segmentations = list(native = bad)), "binding")
    shuffled <- wp02_partition(molecule_id = c("m3", "m1", "m2"),
        status = c("ambiguous", "assigned", "unassigned"),
        cell_id = c(NA_character_, "c1", NA_character_))
    x <- wp03_set(segmentations = list(native = shuffled))
    expect_identical(x@segmentations$native@assignment_fingerprint,
                     wp02_partition()@assignment_fingerprint)
    expect_silent(.checkSetIntegrity(x))
})

test_that("public set data, nested run metadata and structural methods are lossless", {
    m <- wp02_molecules()
    m$coordinate <- c(0, 1, 2)
    x <- wp03_set(molecules = m,
        coordinate_info = S4Vectors::SimpleList(units = "micrometres"),
        provenance = S4Vectors::SimpleList(source = "fixture"),
        metadata = S4Vectors::SimpleList(note = "retained"))
    expect_identical(moleculeData(x), m)
    info <- segmentationInfo(x)
    expect_s4_class(info, "DataFrame")
    expect_identical(colnames(info),
                     c("name", "n_cells", "n_molecules", "run_metadata"))
    expect_identical(info$name, segmentationNames(x))
    expect_identical(info$n_cells, c(3L, 2L))
    expect_identical(info$n_molecules, c(3L, 3L))
    for (i in seq_len(nrow(info))) {
        expect_identical(info$run_metadata[[i]],
                         runMetadata(x@segmentations[[i]]))
    }
    expect_identical(provenance(x)$source, "fixture")
    expect_identical(provenance(x)$software_version, NA_character_)
    expect_identical(x@coordinate_info$units, "micrometres")
    expect_identical(x@metadata$note, "retained")
    s <- summary(x)
    expect_s4_class(s, "SimpleList")
    expect_identical(names(s), c("specimen_id", "n_molecules", "n_features",
        "n_segmentations", "segmentations"))
    expect_identical(s$n_molecules, 3L)
    expect_identical(s$n_features, 2L)
    expect_identical(s$n_segmentations, 2L)
    expect_identical(s$segmentations, info)
    expect_output(methods::show(x), "SegmentationSet")
    expect_output(methods::show(x), "native, other")
    expect_error(runMetadata(x))
    expect_error(cellCensus(x))
    m$molecule_id[1L] <- "changed"
    expect_identical(moleculeData(x)$molecule_id[1L], "m1")
    returned <- moleculeData(x)
    returned$feature_id[1L] <- "changed"
    returned$molecule_id[1L] <- "changed"
    f <- featureData(x)
    f$feature_id[1L] <- "changed"
    p <- provenance(x)
    p$source <- "changed"
    info$run_metadata[[1L]]$nested$parameters$a <- 99L
    expect_identical(moleculeData(x)$feature_id[1L], "g1")
    expect_identical(moleculeData(x)$molecule_id[1L], "m1")
    expect_identical(featureData(x)$feature_id[1L], "g1")
    expect_identical(provenance(x)$source, "fixture")
    expect_identical(segmentationInfo(x)$run_metadata[[1L]]$nested$parameters$a,
                     1L)
    expect_silent(.checkSetIntegrity(x))
})

test_that("cheap validity rejects malformed structure while integrity checks hashes", {
    x <- wp03_set()
    mutations <- list(
        function(z) { z@specimen_id <- character(); z },
        function(z) { z@segmentations <- S4Vectors::SimpleList(); z },
        function(z) { names(z@segmentations) <- c("native", "native"); z },
        function(z) { names(z@segmentations)[1L] <- ""; z },
        function(z) { names(z@segmentations)[1L] <- NA_character_; z },
        function(z) { names(z@segmentations)[1L] <- "alias"; z },
        function(z) { z@segmentations[[1L]] <- 1L; z },
        function(z) { z@segmentations[[1L]]@status <- 1L; z },
        function(z) { z@segmentations[[1L]]@cell_index[1L] <- 99L; z },
        function(z) { z@molecules$molecule_id[2L] <- "m1"; z },
        function(z) { z@molecules$molecule_id <- as.double(1:3); z },
        function(z) { z@molecules$feature_index[1L] <- 99L; z },
        function(z) { z@molecules$feature_index[1L] <- NA_integer_; z },
        function(z) { z@molecules$feature_index <- NULL; z },
        function(z) { z@molecules$feature_index <- c(2L, 1L, 2L); z },
        function(z) { z@molecules$molecule_index <- 3:1; z },
        function(z) { z@features$feature_id[2L] <- "g1"; z },
        function(z) { z@fingerprint_version <- "2"; z },
        function(z) { z@fingerprint_algorithm <- "other"; z },
        function(z) { z@segmentations[[1L]]@fingerprint_version <- "2"; z },
        function(z) { z@segmentations[[1L]]@assignment_fingerprint <- ""; z },
        function(z) { z@schema_version <- ""; z },
        function(z) { z@universe_fingerprint <- ""; z },
        function(z) { z@universe_fingerprint <- character(); z },
        function(z) {
            z@segmentations[[1L]]@universe_fingerprint <-
                paste(rep("0", 64L), collapse = "")
            z
        })
    for (mutate in mutations) {
        expect_error(methods::validObject(mutate(x)))
    }
    bad <- x
    bad@molecules$molecule_id[1L] <- "changed"
    expect_true(methods::validObject(bad))
    expect_error(.checkSetIntegrity(bad), "Universe fingerprint")
    bad <- x
    bad@segmentations[[1L]]@cell_index[1L] <- 2L
    expect_true(methods::validObject(bad))
    expect_error(.checkSetIntegrity(bad), "assignment fingerprint")
})

test_that("historical empty populations remain valid in nonempty sets", {
    m <- S4Vectors::DataFrame(molecule_id = integer(), feature_id = character())
    p <- SegmentationPartition(m, "empty",
        S4Vectors::DataFrame(cell_id = integer()), integer(), character())
    x <- SegmentationSet(m, list(empty = p), "empty-specimen")
    expect_true(methods::validObject(x))
    expect_identical(moleculeData(x), m)
    expect_identical(nrow(featureData(x)), 0L)
    expect_identical(segmentationInfo(x)$n_cells, 0L)
    expect_silent(.checkSetIntegrity(x))
    m <- S4Vectors::DataFrame(molecule_id = 1:2, feature_id = c("g", "g"))
    p <- SegmentationPartition(m, "empty-census",
        S4Vectors::DataFrame(cell_id = integer()), 1:2,
        c("unassigned", "ambiguous"))
    x <- SegmentationSet(m, list("empty-census" = p), "specimen")
    expect_true(methods::validObject(x))
    expect_identical(segmentationInfo(x)$n_molecules, 2L)
})

test_that("set RDS preserves all supported identity families and metadata", {
    for (kind in c("integer", "character", "LLint", "integer64")) {
        if (kind == "integer64" && !requireNamespace("bit64", quietly = TRUE)) {
            next
        }
        ids <- id_lab_fixture(kind)$x
        n <- length(ids)
        m <- S4Vectors::DataFrame(molecule_id = ids,
            feature_id = rep(c("h", "g"), length.out = n))
        p <- SegmentationPartition(m, "run1",
            S4Vectors::DataFrame(cell_id = ids), ids[n:1],
            rep("assigned", n), ids[n:1],
            run_metadata = S4Vectors::SimpleList(seed = 3L))
        q <- p
        q@name <- "run2"
        x <- SegmentationSet(m, list(run2 = q, run1 = p), "specimen",
            coordinate_info = S4Vectors::SimpleList(units = "unknown"),
            provenance = S4Vectors::SimpleList(source = "fixture"),
            metadata = S4Vectors::SimpleList(nested = list(value = ids)))
        path <- tempfile()
        saveRDS(x, path)
        restored <- readRDS(path)
        unlink(path)
        expect_true(methods::validObject(restored))
        expect_id_lab_exact(moleculeData(restored)$molecule_id, ids)
        expect_id_lab_exact(cellCensus(restored@segmentations$run1)$cell_id, ids)
        expect_identical(moleculeData(restored)$feature_id, m$feature_id)
        expect_identical(featureData(restored), featureData(x))
        expect_identical(segmentationNames(restored), c("run2", "run1"))
        expect_identical(restored@universe_fingerprint, x@universe_fingerprint)
        expect_identical(restored@segmentations, x@segmentations)
        expect_identical(restored@coordinate_info, x@coordinate_info)
        expect_identical(restored@provenance, x@provenance)
        expect_identical(restored@metadata, x@metadata)
        expect_silent(.checkSetIntegrity(restored))
    }
})

test_that("construction fingerprints the shared universe once", {
    m <- wp02_molecules()
    p <- wp02_partition()
    q <- wp02_partition(name = "other")
    counter <- new.env(parent = emptyenv())
    counter$n <- 0L
    ns <- asNamespace("SegStability")
    trace(".universeFingerprint", where = ns, print = FALSE,
        tracer = function() counter$n <- counter$n + 1L)
    on.exit(untrace(".universeFingerprint", where = ns), add = TRUE)
    x <- SegmentationSet(m, list(native = p, other = q), "specimen")
    expect_identical(counter$n, 1L)
    expect_true(methods::validObject(x))
    expect_identical(counter$n, 1L)
    .checkSetIntegrity(x)
    expect_identical(counter$n, 2L)
    p@cell_index[1L] <- 2L
    expect_identical(x@segmentations$native@cell_index,
                     c(1L, NA_integer_, NA_integer_))
})
