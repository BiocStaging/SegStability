# Assisted-by: OpenAI Codex
# Development-time replay only. Never sourced by package loading, build, check,
# examples, tests, or vignette rendering. No network access or downloads.
#
# 1. Obtain the exact public members listed in prepare_ovrlpy_validation.py.
# 2. Run: python3 prepare_ovrlpy_validation.py INPUT_DIR WORK_DIR
# 3. Install SegStability 0.99.0, then run:
#      Rscript reproduce_ovrlpy_validation.R WORK_DIR
# Requires the optional bit64 package for original UInt64 transcript IDs.
# Work files stay outside the package; no Xenium data are bundled.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[[1L]])) {
    stop("Usage: Rscript reproduce_ovrlpy_validation.R WORK_DIR",
         call. = FALSE)
}
work_dir <- normalizePath(args[[1L]], mustWork = TRUE)
if (!requireNamespace("bit64", quietly = TRUE)) {
    stop("Install the optional bit64 package for exact Xenium UInt64 IDs",
         call. = FALSE)
}
if (!requireNamespace("SegStability", quietly = TRUE)) {
    stop("Install SegStability before running this development-time replay",
         call. = FALSE)
}
record_file <- system.file("extdata", "validation",
    "ovrlpy_validation_summary.csv", package = "SegStability")
if (!nzchar(record_file)) {
    stop("Installed SegStability lacks the ovrlpy validation record",
         call. = FALSE)
}
record <- utils::read.csv(record_file, stringsAsFactors = FALSE,
                          colClasses = "character")
if (anyDuplicated(record$key) || !identical(names(record),
                                             c("key", "value", "evidence"))) {
    stop("Invalid or duplicate validation-record keys", call. = FALSE)
}
expected <- stats::setNames(record$value, record$key)
check <- function(key, observed) {
    target <- expected[[key]]
    if (is.null(target) || length(observed) != 1L || is.na(observed) ||
        !identical(as.character(observed), target)) {
        stop("Validation mismatch for ", key, ": observed ",
             paste(observed, collapse = ", "), "; expected ", target,
             call. = FALSE)
    }
}
check("segstability_version", as.character(utils::packageVersion("SegStability")))
manifest_file <- file.path(work_dir, "input_verification.csv")
if (!file.exists(manifest_file)) {
    stop("Missing input_verification.csv; run the Python preparation script first",
         call. = FALSE)
}
manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE,
                            colClasses = "character")
if (!identical(names(manifest), c("key", "value")) ||
    anyDuplicated(manifest$key)) {
    stop("Invalid preparation manifest", call. = FALSE)
}
for (i in seq_len(nrow(manifest))) check(manifest$key[[i]], manifest$value[[i]])
read_tab <- function(name, classes) {
    path <- file.path(work_dir, name)
    if (!file.exists(path)) stop("Missing prepared ROI table: ", path,
                                 call. = FALSE)
    utils::read.delim(path, colClasses = classes, check.names = FALSE,
                      stringsAsFactors = FALSE)
}
source <- read_tab("roi_vendor.tsv",
    c(molecule_id = "character", feature_id = "character",
      x_location = "numeric", y_location = "numeric", z_location = "numeric",
      qv = "numeric", vendor_cell_id = "integer"))
baysor <- read_tab("roi_baysor.tsv",
    c(molecule_id = "character", status = "character", cell_id = "character"))
source_cells <- read_tab("roi_vendor_cells.tsv",
    c(cell_id = "integer", x_centroid = "numeric", y_centroid = "numeric"))
baysor_cells <- read_tab("roi_baysor_cells.tsv",
    c(cell_id = "character", x = "numeric", y = "numeric"))
check("roi_universe_size", nrow(source))
check("roi_universe_size", nrow(baysor))
check("vendor_census_cells", nrow(source_cells))
check("baysor_census_cells", nrow(baysor_cells))
if (!identical(source$molecule_id, baysor$molecule_id) ||
    anyDuplicated(source$molecule_id) ||
    any(!baysor$status %in% c("assigned", "unassigned")) ||
    anyNA(source$feature_id) || any(!nzchar(source$feature_id))) {
    stop("Prepared ROI IDs, features, or statuses differ", call. = FALSE)
}
if (anyNA(baysor$cell_id[baysor$status == "unassigned"]) ||
    any(baysor$cell_id[baysor$status == "unassigned"] != "")) {
    stop("Unassigned Baysor molecules must have blank cell IDs", call. = FALSE)
}
baysor$cell_id[baysor$status == "unassigned"] <- NA_character_
ids <- bit64::as.integer64(source$molecule_id)
if (!identical(as.character(ids), source$molecule_id)) {
    stop("Original molecule IDs failed integer64 decimal round-trip",
         call. = FALSE)
}
molecules <- data.frame(molecule_id = ids, feature_id = source$feature_id,
    x_location = source$x_location, y_location = source$y_location,
    z_location = source$z_location, qv = source$qv)
vendor_cells <- data.frame(cell_id = source_cells$cell_id,
    x_centroid = source_cells$x_centroid,
    y_centroid = source_cells$y_centroid)
baysor_cells <- data.frame(cell_id = baysor_cells$cell_id,
    x = baysor_cells$x, y = baysor_cells$y)
vendor_id <- source$vendor_cell_id
vendor_id[vendor_id == -1L] <- NA_integer_
vendor <- SegStability::SegmentationPartition(molecules, "vendor",
    vendor_cells, ids, ifelse(is.na(vendor_id), "unassigned", "assigned"),
    vendor_id, run_metadata = S4Vectors::SimpleList(
        software = "Xenium", version = "1.0.2"))
baysor_partition <- SegStability::SegmentationPartition(molecules, "baysor",
    baysor_cells, bit64::as.integer64(baysor$molecule_id), baysor$status,
    baysor$cell_id, run_metadata = S4Vectors::SimpleList(
        software = "Baysor", version = "0.7.1", run_id = expected[["baysor_run_id"]]))
set <- SegStability::SegmentationSet(molecules,
    list(vendor = vendor, baysor = baysor_partition),
    specimen_id = paste0(expected[["experiment"]], ":replicate1:tile9x9"))
if (!isTRUE(methods::validObject(set, test = TRUE))) {
    stop("Constructed SegmentationSet is invalid", call. = FALSE)
}
comparison <- SegStability::compareSegmentations(set, "vendor", "baysor")
if (!isTRUE(methods::validObject(comparison, test = TRUE))) {
    stop("Constructed SegmentationComparison is invalid", call. = FALSE)
}
edge <- SegStability::correspondenceEdges(comparison)
left <- SegStability::cellSummary(comparison, "vendor")
right <- SegStability::cellSummary(comparison, "baysor")
components <- SegStability::componentSummary(comparison)
status <- SegStability::assignmentSummary(comparison)
check("vendor_assigned", sum(left$assigned_mass))
check("baysor_assigned", sum(right$assigned_mass))
check("vendor_unassigned", sum(source$vendor_cell_id == -1L))
check("baysor_unassigned", sum(baysor$status == "unassigned"))
check("vendor_zero_roi_cells", sum(left$assigned_mass == 0))
check("baysor_zero_roi_cells", sum(right$assigned_mass == 0))
check("correspondence_edges", length(edge))
check("shared_mass", sum(S4Vectors::mcols(edge)$weight))
check("components", nrow(components))
check("isolate_components", sum(components$n_segmentation_1 +
    components$n_segmentation_2 == 1L))
check("many_to_many_components", sum(components$component_class ==
    "many_to_many"))
check("vendor_unique_dominance_cells", sum(left$dominance_count == 1L))
check("vendor_tied_dominance_cells", sum(left$dominance_count > 1L))
check("baysor_unique_dominance_cells", sum(right$dominance_count == 1L))
check("baysor_tied_dominance_cells", sum(right$dominance_count > 1L))
check("reciprocal_dominant_edges",
      sum(S4Vectors::mcols(edge)$reciprocal_dominant))
for (a in c("assigned", "unassigned", "ambiguous")) {
    for (b in c("assigned", "unassigned", "ambiguous")) {
        key <- paste("status", a, b, sep = "_")
        if (key %in% names(expected)) {
            row <- status$status_1 == a & status$status_2 == b
            check(key, status$n_molecules[row])
        }
    }
}
thresholds <- c(0, 0.01, 0.025, 0.05, 0.1, 0.2)
labels <- c("0", "0.01", "0.025", "0.05", "0.1", "0.2")
support <- SegStability::supportSensitivity(comparison,
    thresholds = thresholds, direction = NULL, states = NULL)
for (side in c("vendor", "baysor")) {
    for (i in seq_along(thresholds)) {
        row <- which(support$direction == side &
                     support$threshold == thresholds[[i]])
        if (length(row) != 1L) stop("Missing directional support row",
                                   call. = FALSE)
        key <- paste0("support_", side, "_", labels[[i]])
        check(paste0(key, "_retained_edges"), support$retained_edges[row])
        check(paste0(key, "_retained_mass"), support$retained_mass[row])
        check(if (side == "vendor") "vendor_assigned" else "baysor_assigned",
              support$assigned_mass[row])
        check("shared_mass", support$original_matched_mass[row])
        if (support$retained_mass[row] +
            support$excluded_matched_mass[row] != support$original_matched_mass[row]) {
            stop("Support mass accounting differs", call. = FALSE)
        }
        detail <- support$cells[[row]]
        if (!identical(detail$vendor$assigned_mass, left$assigned_mass) ||
            !identical(detail$baysor$assigned_mass, right$assigned_mass) ||
            !identical(detail$vendor$component_id, left$component_id) ||
            !identical(detail$baysor$component_id, right$component_id) ||
            !identical(support$components[[row]]$component_class,
                       components$component_class)) {
            stop("Support filtering altered original denominators/components",
                 call. = FALSE)
        }
    }
}
cat("Public-source replay agrees with the bundled record: ",
    nrow(source), " molecules, ", length(edge), " correspondence edges, ",
    nrow(support), " directional support rows.\n", sep = "")
