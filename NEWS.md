# SegStability 0.99.0

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
