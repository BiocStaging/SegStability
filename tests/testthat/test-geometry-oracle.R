test_that("validated historical raw polygons reproduce saved overlap metrics", {
    skip_if_not_installed("sf")
    cmp <- wp04_compare(matrix(5L, 1, 1), ids1 = "102",
                         ids2 = "CR3b816afa8-388")
    shapes <- wp09_oracle_shapes()
    g <- wp09_geometry(cmp, shapes[1L], shapes[2L])
    out <- wp09_run(cmp, g)
    # Saved row of fix_phase7/geometry_comparison/baysor_edge_geometry.csv.gz.
    expect_equal(out$pairs$intersection_area, 5.42672533921602, tolerance = 1e-10)
    expect_equal(out$pairs$iou, .010479030392815296, tolerance = 1e-10)
    expect_equal(out$pairs$dice, .020740718169562928, tolerance = 1e-10)
    expect_equal(out$pairs$overlap_fraction_1, .013171478556623649, tolerance = 1e-10)
    expect_equal(out$pairs$overlap_fraction_2, .048763704541261305, tolerance = 1e-10)
    expect_identical(out$pairs$cell_id_1, "102")
    expect_identical(out$pairs$cell_id_2, "CR3b816afa8-388")
})
