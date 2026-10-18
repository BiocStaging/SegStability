# Assisted-by: OpenAI Codex
## Regenerate data/segmentationExample.rda with Rscript and an output path.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
    stop("Supply the path to segmentationExample.rda")
}

ids <- sprintf("m%02d", seq_len(24L))
first_cell <- c(rep("a1", 4L), rep("a2", 6L), rep("a3", 3L),
                rep("a4", 3L), rep("a5", 4L), "a1", NA, NA, NA)
second_cell <- c(rep("b1", 4L), rep("b2", 3L), rep("b3", 3L),
                 rep("b4", 6L), rep("b5", 2L), rep("b6", 2L),
                 NA, "b1", "b2", NA)
first_status <- c(rep("assigned", 21L), "unassigned", "ambiguous",
                  "unassigned")
second_status <- c(rep("assigned", 20L), "unassigned", "assigned",
                   "assigned", "unassigned")

segmentationExample <- list(
    molecules = data.frame(
        molecule_id = ids,
        feature_id = rep(c("gA", "gB", "gC"), length.out = length(ids)),
        x = as.numeric(rep(seq_len(6L), each = 4L)),
        y = as.numeric(rep(seq_len(4L), times = 6L))),
    cells = list(
        first = data.frame(cell_id = c(paste0("a", 1:5), "a0")),
        second = data.frame(cell_id = c(paste0("b", 1:6), "b0"))),
    assignments = list(
        first = data.frame(molecule_id = ids, status = first_status,
                           cell_id = first_cell),
        second = data.frame(molecule_id = ids, status = second_status,
                            cell_id = second_cell)),
    ambiguity = list(first = list(
        evidence = data.frame(molecule_id = "m23", candidate_count = 2L,
                              candidate_ids_known = TRUE),
        candidates = data.frame(molecule_id = c("m23", "m23"),
                                cell_id = c("a1", "a2"))),
        second = NULL),
    states = list(
        first = c(a1 = "s1", a2 = "s2", a3 = "s3", a4 = "s3",
                  a5 = "s4", a0 = NA_character_),
        second = c(b1 = "s1", b2 = "s2", b3 = "s5", b4 = "s3",
                   b5 = "s4", b6 = NA_character_, b0 = NA_character_)))

save(segmentationExample, file = args[[1L]], compress = "xz", version = 2L)
