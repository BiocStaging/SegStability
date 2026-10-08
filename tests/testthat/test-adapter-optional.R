test_that("core namespace loads without loading optional adapter namespaces", {
    expect_true(isNamespaceLoaded("SegStability"))
    if (!requireNamespace("MoleculeExperiment", quietly = TRUE)) {
        expect_error(fromMoleculeExperiment(NULL), "Optional package MoleculeExperiment")
    }
    if (!requireNamespace("SpatialFeatureExperiment", quietly = TRUE)) {
        expect_error(geometryFromSpatialFeatureExperiment(NULL),
            "Optional package SpatialFeatureExperiment")
    }
    expect_s4_class(wp11_manual(wp11_tables()), "SegmentationSet")
})
