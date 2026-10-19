# Assisted-by: OpenAI Codex

test_that("bundled ovrlpy summary is parseable and internally consistent", {
    path <- system.file("extdata", "validation",
        "ovrlpy_validation_summary.csv", package = "SegStability")
    expect_true(nzchar(path))
    record <- read.csv(path, colClasses = "character")
    expect_named(record, c("key", "value", "evidence"))
    expect_false(anyDuplicated(record$key) > 0L)
    expect_false(anyNA(record$key) || anyNA(record$value))
    values <- stats::setNames(record$value, record$key)
    count <- function(key) as.numeric(values[[key]])

    expect_equal(count("status_assigned_assigned") +
        count("status_assigned_unassigned") +
        count("status_unassigned_assigned") +
        count("status_unassigned_unassigned"), count("roi_universe_size"))
    expect_equal(count("status_assigned_assigned") +
        count("status_assigned_unassigned"), count("vendor_assigned"))
    expect_equal(count("status_assigned_assigned") +
        count("status_unassigned_assigned"), count("baysor_assigned"))
    expect_equal(count("vendor_assigned") + count("vendor_unassigned") +
        count("vendor_ambiguous"), count("roi_universe_size"))
    expect_equal(count("baysor_assigned") + count("baysor_unassigned") +
        count("baysor_ambiguous"), count("roi_universe_size"))
    expect_equal(count("shared_mass"), count("status_assigned_assigned"))

    for (side in c("vendor", "baysor")) {
        levels <- c("0", "0.01", "0.025", "0.05", "0.1", "0.2")
        edges <- vapply(levels, function(level)
            count(paste0("support_", side, "_", level,
                         "_retained_edges")), numeric(1))
        mass <- vapply(levels, function(level)
            count(paste0("support_", side, "_", level,
                         "_retained_mass")), numeric(1))
        expect_equal(edges[[1L]], count("correspondence_edges"))
        expect_equal(mass[[1L]], count("shared_mass"))
        expect_true(all(diff(edges) <= 0))
        expect_true(all(diff(mass) <= 0))
        expect_true(all(mass <= count("shared_mass")))
    }
    hashes <- values[grepl("^sha256_", names(values))]
    expect_true(all(grepl("^[0-9a-f]{64}$", hashes)))
})
