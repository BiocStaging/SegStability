validation_row <- function(x, check) x[x$check == check, , drop = FALSE]

validation_defect <- function(x) any(x$status == "FAIL" & x$severity == "error")
