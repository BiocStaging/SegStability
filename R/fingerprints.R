## Version 1: length-framed fields with explicit family and missingness tags.
## Canonical text uses UTF-8 without normalization; byte-marked IDs stay bytes.
## Ancillary metadata/provenance are deliberately outside identity binding.
.frameValues <- function(x) {
    if (!length(x)) {
        return("")
    }
    missing <- is.na(x)
    x[missing] <- ""
    byte <- Encoding(x) == "bytes"
    x[!byte] <- enc2utf8(x[!byte])
    paste0(ifelse(missing, "N", ifelse(byte, "B", "V")),
           nchar(x, type = "bytes"), ":", x, collapse = "")
}

.fingerprint <- function(domain, fields) {
    records <- vapply(seq_along(fields), function(i) {
        x <- fields[[i]]
        paste0(.frameValues(names(fields)[i]), length(x), ":",
               .frameValues(unname(as.character(x))))
    }, character(1))
    bytes <- charToRaw(paste0("SegStability:1:", .frameValues(domain),
                             paste0(records, collapse = "")))
    S4Vectors::SimpleList(fingerprint_version = "1",
        fingerprint_algorithm = "sha256",
        fingerprint_value = digest::digest(bytes, algo = "sha256",
                                            serialize = FALSE))
}

.universeFingerprint <- function(molecules, features) {
    .fingerprint("universe", list(
        molecule_family = .idFamily(molecules$molecule_id),
        molecule_id = .idKeys(molecules$molecule_id),
        feature_id = features$feature_id,
        feature_index = molecules$feature_index))
}

.assignmentFingerprint <- function(universe_fingerprint, cells, status,
                                    cell_index, ambiguity) {
    .fingerprint("assignment", list(
        universe_fingerprint = universe_fingerprint,
        cell_family = .idFamily(cells$cell_id),
        cell_id = .idKeys(cells$cell_id), status = status,
        cell_index = cell_index,
        ambiguous_molecule = ambiguity$evidence$molecule_index,
        candidate_count = ambiguity$evidence$candidate_count,
        candidate_ids_known = ambiguity$evidence$candidate_ids_known,
        candidate_molecule = S4Vectors::from(ambiguity$candidates),
        candidate_cell = S4Vectors::to(ambiguity$candidates)))
}
