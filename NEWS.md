# SegStability 0.99.0

* Represent alternative segmentations of one molecule universe with exact IDs,
  complete cell censuses, and explicit assigned, unassigned, and ambiguous
  statuses.
* Describe shared-molecule correspondence, categorical cell-state stability,
  directional support sensitivity, optional topology and geometry, and
  agreement among segmentation runs.
* Support explicit MoleculeExperiment molecule import and keyed
  SpatialFeatureExperiment geometry sidecars when those packages are installed.
* Retrieve a stored partition with `set[["name"]]`, read its keyed assignments
  with `assignmentData(set, "name")`, and inspect stored coordinate declarations
  with `summary(set)$coordinate_info` or
  `summary(comparison)$coordinate_info`.
