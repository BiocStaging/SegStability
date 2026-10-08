.ingestionTable <- function(x, label) {
    if (methods::is(x, "DataFrame")) return(x)
    if (!inherits(x, "data.frame") || inherits(x, "sf")) {
        stop(label, " must be a data.frame or DataFrame", call. = FALSE)
    }
    if (anyDuplicated(names(x))) {
        stop(label, " has duplicate column names", call. = FALSE)
    }
    out <- S4Vectors::DataFrame(x, check.names = FALSE)
    for (j in seq_along(x)) {
        if (!identical(out[[j]], x[[j]])) {
            stop(label, ": column representation changed during ingestion",
                 call. = FALSE)
        }
    }
    out
}

.adapterPackage <- function(package, object) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Optional package ", package, " is required for this adapter",
             call. = FALSE)
    }
    if (!isS4(object) || !methods::is(object, package)) {
        stop("Expected a ", package, " source object", call. = FALSE)
    }
}

.adapterScalar <- function(x, label) {
    if (!.scalarString(x)) {
        stop(label, " must be a nonempty character scalar", call. = FALSE)
    }
    invisible(x)
}

.adapterValidation <- function(report) {
    bad <- report$status == "FAIL" & report$severity == "error"
    if (any(bad)) {
        stop(paste(report$message[bad], collapse = "; "), call. = FALSE)
    }
    invisible(report)
}

.adapterSpecimen <- function(coordinate_info, specimen_id) {
    problems <- .coordinateListProblems(coordinate_info)
    if (!isTRUE(problems)) stop(problems, call. = FALSE)
    declared <- coordinate_info$specimen_id
    if (!is.null(declared) && !is.na(declared) &&
        !identical(declared, specimen_id)) {
        stop("coordinate_info specimen_id conflicts with specimen_id",
             call. = FALSE)
    }
}

.adapterSourceSpecimen <- function(molecules, specimen_id) {
    declared <- molecules$specimen_id
    if (is.null(declared)) return(invisible(NULL))
    if (!is.character(declared) || is.object(declared) ||
        any(declared[!is.na(declared)] != specimen_id)) {
        stop("Source specimen_id conflicts with explicit specimen_id",
             call. = FALSE)
    }
    invisible(NULL)
}

.meCoordinates <- function(table, coordinate_cols) {
    if (is.null(coordinate_cols)) return(table)
    axes <- names(coordinate_cols)
    if (!is.character(coordinate_cols) || is.null(axes) ||
        anyNA(coordinate_cols) || any(!nzchar(coordinate_cols)) ||
        anyDuplicated(coordinate_cols) ||
        !identical(axes, c("x", "y")) &&
        !identical(axes, c("x", "y", "z"))) {
        stop("coordinate_cols must map x/y or x/y/z to distinct columns",
             call. = FALSE)
    }
    if (!all(coordinate_cols %in% names(table))) {
        stop("Selected coordinate columns are missing", call. = FALSE)
    }
    for (j in seq_along(axes)) {
        source <- coordinate_cols[j]
        axis <- axes[j]
        if (axis %in% names(table) && source != axis &&
            !identical(table[[axis]], table[[source]])) {
            stop("Conflicting source coordinate columns", call. = FALSE)
        }
    }
    values <- lapply(coordinate_cols, function(column) table[[column]])
    for (source in setdiff(coordinate_cols, axes)) table[[source]] <- NULL
    for (j in seq_along(axes)) table[[axes[j]]] <- values[[j]]
    table
}

.meMoleculeTable <- function(leaf, feature_id, sample_id, molecule_id_col,
                             coordinate_cols) {
    t <- .ingestionTable(leaf, "source molecule table")
    if (anyDuplicated(names(t)) || !molecule_id_col %in% names(t)) {
        stop("Original molecule ID column is missing or nonunique",
             call. = FALSE)
    }
    ids <- .validateIds(t[[molecule_id_col]], "source molecule ID")
    feature <- rep(feature_id, nrow(t))
    if ("feature_id" %in% names(t) &&
        !identical(.validateFeatures(t$feature_id), feature)) {
        stop("Source feature_id conflicts with retained feature key",
             call. = FALSE)
    }
    if (molecule_id_col != "molecule_id") {
        if ("molecule_id" %in% names(t)) {
            stop("Conflicting molecule_id columns", call. = FALSE)
        }
        names(t)[names(t) == molecule_id_col] <- "molecule_id"
    }
    t$feature_id <- feature
    t <- .meCoordinates(t, coordinate_cols)
    if ("sample_id" %in% names(t) &&
        !identical(t$sample_id, rep(sample_id, nrow(t)))) {
        stop("Source sample_id conflicts with selected sample",
             call. = FALSE)
    }
    t
}

.meMolecules <- function(object, assay, sample_id, molecule_id_col,
                         coordinate_cols) {
    nested <- MoleculeExperiment::molecules(object, assayName = assay,
                                            flatten = FALSE)
    samples <- nested[[assay]]
    if (is.null(names(samples)) || anyDuplicated(names(samples)) ||
        !sample_id %in% names(samples)) {
        stop("Explicit sample_id must select a unique source sample",
             call. = FALSE)
    }
    leaves <- samples[[sample_id]]
    if (!is.list(leaves) || !length(leaves) || is.null(names(leaves))) {
        stop("Source sample has no retained molecule tables", call. = FALSE)
    }
    .validateFeatures(names(leaves), unique = TRUE)
    tables <- lapply(seq_along(leaves), function(j)
        .meMoleculeTable(leaves[[j]], names(leaves)[j], sample_id,
                         molecule_id_col, coordinate_cols))
    .meBindTables(tables)
}

.meBindTables <- function(tables) {
    keys <- vector("list", length(tables))
    family <- NULL
    for (j in seq_along(tables)) {
        t <- tables[[j]]
        ids <- t$molecule_id
        kind <- .idFamily(ids)
        if (!is.null(family) && !identical(kind, family)) {
            stop("Source molecule ID families must agree", call. = FALSE)
        }
        family <- kind
        keys[[j]] <- .idKeys(ids)
        if (j > 1L) {
            if (!identical(names(t), names(tables[[1L]])) ||
                any(!vapply(seq_along(t), function(k)
                    identical(class(t[[k]]), class(tables[[1L]][[k]])),
                    logical(1)))) {
                stop("Incompatible source molecule table schemas",
                     call. = FALSE)
            }
        }
    }
    out <- do.call(rbind, tables)
    .validateIds(out$molecule_id, "source molecule ID")
    if (!identical(.idFamily(out$molecule_id), family) ||
        !identical(.idKeys(out$molecule_id), unlist(keys, use.names = FALSE))) {
        stop("Molecule IDs changed while combining source tables",
             call. = FALSE)
    }
    canonical <- intersect(c("molecule_id", "feature_id", "x", "y", "z"),
                           names(out))
    out <- out[, c(canonical, setdiff(names(out), canonical)), drop = FALSE]
    rownames(out) <- NULL
    out
}

.meUniverse <- function(source, universe) {
    if (is.null(universe)) return(source)
    universe <- .ingestionTable(universe, "universe")
    .prepareUniverse(universe)
    order <- .alignMolecules(source$molecule_id, universe$molecule_id)
    source <- source[order, , drop = FALSE]
    if (!identical(source$feature_id, universe$feature_id)) {
        stop("Source feature identities conflict with universe", call. = FALSE)
    }
    for (column in intersect(names(source), names(universe))) {
        if (!identical(source[[column]], universe[[column]])) {
            label <- if (column %in% c("x", "y", "z")) "coordinate" else column
            stop("Source ", label, " conflicts with universe", call. = FALSE)
        }
    }
    for (column in setdiff(names(source), names(universe))) {
        universe[[column]] <- source[[column]]
    }
    universe
}

.sfeSourceGeometry <- function(object, type, sample_id, geometry_kind) {
    .adapterScalar(type, "type")
    .adapterScalar(sample_id, "sample_id")
    if (!identical(geometry_kind, "col") &&
        !identical(geometry_kind, "annot")) {
        stop("geometry_kind must be col or annot", call. = FALSE)
    }
    geometry <- if (geometry_kind == "col") {
        if (is.null(colnames(object))) {
            stop("Column geometry requires explicit source column keys",
                 call. = FALSE)
        }
        SpatialFeatureExperiment::colGeometry(object, type = type,
            sample_id = sample_id, withDimnames = TRUE)
    } else {
        SpatialFeatureExperiment::annotGeometry(object, type = type,
                                                 sample_id = sample_id)
    }
    if (!inherits(geometry, "sf") || .row_names_info(geometry, 1L) < 0L) {
        stop("Geometry requires explicitly retained source row keys",
             call. = FALSE)
    }
    keys <- rownames(geometry)
    .validateIds(keys, "geometry source keys")
    geometry
}

.sfeMappedGeometry <- function(geometry, cell_map) {
    keys <- rownames(geometry)
    cell_map <- .ingestionTable(cell_map, "cell_map")
    if (anyDuplicated(names(cell_map)) ||
        !all(c("geometry_id", "cell_id") %in% names(cell_map))) {
        stop("cell_map requires geometry_id and cell_id", call. = FALSE)
    }
    .validateFeatures(cell_map$geometry_id, unique = TRUE)
    .validateIds(cell_map$cell_id, "mapped cell ID")
    index <- .alignMolecules(cell_map$geometry_id, keys)
    ids <- cell_map$cell_id[index]
    if ("cell_id" %in% names(geometry) &&
        !identical(geometry$cell_id, ids)) {
        stop("Source cell_id conflicts with explicit cell_map", call. = FALSE)
    }
    geometry$cell_id <- ids
    geometry
}
