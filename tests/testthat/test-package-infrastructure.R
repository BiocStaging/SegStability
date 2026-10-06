test_that("the skeleton loads without scientific exports", {
    expect_true(isNamespaceLoaded("SegStability"))
    expect_identical(getNamespaceExports("SegStability"), character())
    expect_identical(
        as.character(utils::packageVersion("SegStability")), "0.99.0"
    )
})
