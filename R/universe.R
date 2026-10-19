# Assisted-by: OpenAI Codex
.prepareUniverse <- function(molecules) {
    molecules <- .ingestionTable(molecules, "molecules")
    if (!methods::is(molecules, "DataFrame") ||
        !all(c("molecule_id", "feature_id") %in% colnames(molecules)) ||
        anyDuplicated(colnames(molecules))) {
        stop("molecules must be a DataFrame with molecule_id and feature_id",
             call. = FALSE)
    }
    if (any(c("molecule_index", "feature_index") %in% colnames(molecules))) {
        stop("Supply feature_id; internal index columns are reserved",
             call. = FALSE)
    }
    .checkIndexSize(nrow(molecules))
    .validateIds(molecules$molecule_id, "molecule_id")
    feature <- .validateFeatures(molecules$feature_id)
    features <- S4Vectors::DataFrame(feature_id = unique(feature))
    molecules$feature_id <- NULL
    molecules$molecule_index <- seq_len(nrow(molecules))
    molecules$feature_index <- base::match(feature, features$feature_id)
    fp <- .universeFingerprint(molecules, features)
    S4Vectors::SimpleList(molecules = molecules, features = features,
        fingerprint_version = fp$fingerprint_version,
        fingerprint_algorithm = fp$fingerprint_algorithm,
        fingerprint_value = fp$fingerprint_value)
}

.checkUniverseStructure <- function(m, f) {
    if (!methods::is(m, "DataFrame") || !methods::is(f, "DataFrame") ||
        !all(c("molecule_id", "molecule_index", "feature_index") %in%
             colnames(m)) || !"feature_id" %in% colnames(f) ||
        anyDuplicated(colnames(m)) || anyDuplicated(colnames(f))) {
        stop("Invalid internal universe structure", call. = FALSE)
    }
    .checkIndexSize(nrow(m))
    .validateIds(m$molecule_id, "molecule_id")
    .validateFeatures(f$feature_id, unique = TRUE)
    if (!identical(m$molecule_index, seq_len(nrow(m))) ||
        !is.integer(m$feature_index) || length(m$feature_index) != nrow(m) ||
        anyNA(m$feature_index) || any(m$feature_index < 1L) ||
        any(m$feature_index > nrow(f))) {
        stop("Invalid molecule_index or feature_index", call. = FALSE)
    }
    if (!identical(unique(m$feature_index), seq_len(nrow(f)))) {
        stop("Feature dictionary must follow first-occurrence mapping",
             call. = FALSE)
    }
    invisible(NULL)
}

.validateUniverse <- function(universe) {
    if (!methods::is(universe, "SimpleList")) {
        stop("Invalid internal universe structure", call. = FALSE)
    }
    m <- universe$molecules
    f <- universe$features
    .checkUniverseStructure(m, f)
    if (!identical(universe$fingerprint_version, "1") ||
        !identical(universe$fingerprint_algorithm, "sha256")) {
        stop("Unsupported universe fingerprint version/algorithm",
             call. = FALSE)
    }
    if (!identical(universe$fingerprint_value,
                   .universeFingerprint(m, f)$fingerprint_value)) {
        stop("Universe fingerprint mismatch", call. = FALSE)
    }
    invisible(NULL)
}
