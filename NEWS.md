# SegStability 0.99.0

* Expose exact stored coordinate declarations in set and comparison summaries,
  including partial and unknown declarations without spatial prerequisites.
* Add owning-set `assignmentData()` readout with exact molecule IDs, semantic
  statuses and native typed cell IDs in canonical universe order.
* Add exact named retrieval of stored, universe-bound partitions with `[[`.

* Add conditional MoleculeExperiment ingestion with original exact molecule IDs,
  explicit assignments and full censuses, plus explicitly keyed
  SpatialFeatureExperiment geometry sidecars. Preserve ambiguity, coordinates,
  CRS and supplied metadata without inference or repair. Adapter packages remain
  optional; base data.frame ingestion preserves scientific column classes.
* Initialize package infrastructure for descriptive segmentation-sensitivity
  analysis.
* Add exact molecule-universe alignment, character-only feature dictionaries,
  and versioned internal universe and assignment integrity fingerprints.
* Add SegmentationPartition with explicit cell censuses, three assignment
  statuses, retained ambiguity evidence, structural validity and accessors.
* Add SegmentationSet ownership and exact partition binding, ordered segmentation
  descriptors, molecule/feature accessors and structural show/summary methods.
* Add SegmentationComparison and compareSegmentations() with complete positive
  correspondence, independent denominators, exact dominance ties, unique
  reciprocity, components including isolates and assignment-status accounting.
* Add explicit partition, set, coordinate, geometry and provenance validation
  reports separating consistency, completeness, integrity and unchecked evidence.
  Finite double coordinates are supported; strict mode stops only on FAIL/error.
  Optional sf geometry is inspected without repair or transformation.
* Add stateStability() for external categorical labels, with complete-graph
  weighted agreement, disagreement and unresolved accounting, directional
  equal-cell means and numeric component summaries. Preserve supplied state
  types, independent denominators and distinct-molecule ambiguity union.
  Required state checks must PASS; unavailable source evidence stays UNKNOWN.
* Add supportSensitivity() with independent directional fractional thresholds,
  original denominators and graph semantics, full-census cell/component support
  accounting and optional state estimands. Structural results are independent
  of labels; unavailable source evidence remains nonblocking UNKNOWN.
* Add anchorTopology() with complete-graph unique reciprocal anchors, strict
  rectangular-window bbox eligibility, exact secondary fan classes and explicit
  unavailable/ineligible geometry strata. Preserve full censuses and graph
  semantics without geometry repair or coordinate transformation.
* Add compareGeometry() with polygon overlap metrics for molecule edges and
  positive geometric pairs, full-census eligibility and exact-objective derived
  maximum-IoU matching. Preserve primary graph/components and source geometry.
* Add replicateAgreement() for explicitly selected segmentation runs, reporting
  all pairwise assignment-status contingencies and conditional label-invariant
  ARI. Preserve run metadata and exact combinatorial degeneracy checks without
  matching run-local cell labels or averaging run-specific conclusions.
