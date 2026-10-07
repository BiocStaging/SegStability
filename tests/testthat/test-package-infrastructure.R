test_that("the namespace exposes exactly the WP08 public boundary", {
    expect_true(isNamespaceLoaded("SegStability"))
    exports <- getNamespaceExports("SegStability")
    public <- exports[!startsWith(exports, ".__")]
    expect_identical(sort(public), sort(c("SegmentationPartition", "cellCensus",
        "assignmentStatus", "ambiguityData", "provenance", "runMetadata",
        "show", "SegmentationSet", "segmentationNames", "segmentationInfo",
        "moleculeData", "featureData", "summary", "compareSegmentations",
        "correspondenceEdges", "cellSummary", "componentSummary",
        "assignmentSummary", "validateSegmentation", "validateSegmentationSet",
        "validateCoordinates", "validateGeometry", "validateProvenance",
        "stateStability", "supportSensitivity", "anchorTopology")))
    expect_identical(
        as.character(utils::packageVersion("SegStability")), "0.99.0"
    )
})

test_that("the S4 show generic is explicitly imported for minimal namespace load", {
    imports <- getNamespaceImports("SegStability")
    expect_true("show" %in% names(imports$methods))
})
