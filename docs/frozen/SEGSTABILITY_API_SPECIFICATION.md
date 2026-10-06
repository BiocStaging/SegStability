# SEGSTABILITY_API_SPECIFICATION.md

**Status: FROZEN FOR v0.99 IMPLEMENTATION**  
**Freeze date:** 2026-10-05  
**Architecture decision:** READY TO FREEZE  

This document freezes the public API and internal object model for SegStability v0.99.x. It implements the scientific contract in `SEGSTABILITY_SCIENTIFIC_SPECIFICATION.md` and the final pre-API decision in `SEGSTABILITY_FINAL_PRE_API_AUDIT.md`.

No redesign is permitted during implementation unless a concrete contradiction is demonstrated by implementation, testing, or Bioconductor compatibility.

SegStability is a descriptive segmentation-sensitivity analytical layer. It is not a segmentation algorithm, segmentation-quality truth score, biological-ground-truth system, cell-type caller, consensus-segmentation method, or replacement spatial-data container.

---

## 1. Frozen architectural principles

SegStability owns exactly three S4 classes in v0.99:

1. `SegmentationPartition`
2. `SegmentationSet`
3. `SegmentationComparison`

Existing Bioconductor infrastructure is reused for storage:

- `S4Vectors::DataFrame`
- `S4Vectors::Hits`
- `S4Vectors::SimpleList`
- `S4Vectors::LLint` where validated for lossless 64-bit identifiers

No additional package-owned class is required for provenance, geometry, states, support sensitivity, replicate results, or graphs in v0.99.

Geometry is optional and external to the core S4 model.

The complete weighted molecule-correspondence graph and independent denominators are primary. Geometry, biological states, support thresholds, reciprocal-anchor topology, and replicate analyses are derived layers and must never redefine the primary graph.

---

## 2. Lossless identifier contract

### 2.1 Scientific semantics

Molecule IDs and cell IDs are opaque identifiers. Equal printed values across unrelated segmentation namespaces do not establish identity.

Supported identifier families may include:

- base R `integer`
- base R `character`
- `S4Vectors::LLint`
- `bit64::integer64`
- UUID-like character identifiers

### 2.2 Conditional support for LLint/integer64

`LLint` and `integer64` belong to the lossless-ID contract, but implementation support is conditional on explicit tests demonstrating exact behavior for:

- duplicate detection
- equality/matching
- ordering
- subsetting
- `S4Vectors::DataFrame` storage and extraction
- serialization/deserialization round trips
- joins/reordering used by universe alignment

No fallback through ordinary R double/numeric is permitted.

If either representation cannot satisfy these tests exactly in the supported Bioconductor/R environment, that representation must fail clearly rather than be coerced through double.

### 2.3 Rejected identifiers

Reject:

- base R `double`/`numeric` identifiers
- logical identifiers
- complex/raw identifiers
- missing identifiers
- unsupported factor identifiers unless explicitly and losslessly converted by the user before construction
- arbitrary unsupported classes

SegStability does not guess that a floating-point value represents a safe integer identity.

### 2.4 Internal computational identity

After exact identity verification, internal computation uses compact positional indices:

- `molecule_index = 1..|U|`
- `cell_index = 1..n_cells` within each segmentation

Molecule position is the canonical computational key after universe alignment. Cell indices are segmentation-local.

Original scientific identifiers remain available once, in the universe/census, for audit, joins, export, diagnostics, and provenance.

Original IDs are not duplicated into correspondence edges.

---

## 3. Safe alignment to the immutable molecule universe

A `SegmentationPartition` is always universe-bound.

Its constructor receives the universe and external molecule IDs for the assignment data. Construction must:

1. validate identifier representation;
2. detect duplicate assignment molecule IDs;
3. compare assignment identity against the immutable universe;
4. reject missing universe molecules;
5. reject extra molecules;
6. calculate exact identity-based ordering;
7. reorder assignment information exactly once into universe order;
8. translate assigned cell IDs into segmentation-local `cell_index`;
9. discard duplicated assignment-side molecule-ID storage after binding;
10. store a universe fingerprint and assignment fingerprint.

Length equality is never sufficient evidence of identity alignment.

After successful construction, assignment vectors are position-aligned to the immutable universe and downstream computation may safely be index based.

No separate assignment S4 class is introduced.

---

## 4. Fingerprints and serialization versioning

Every `SegmentationPartition` stores:

- `universe_fingerprint`
- `assignment_fingerprint`
- `fingerprint_version`
- `fingerprint_algorithm`

Every `SegmentationSet` stores:

- `universe_fingerprint`
- `fingerprint_version`
- `fingerprint_algorithm`

Every `SegmentationComparison` stores:

- `universe_fingerprint`
- source `assignment_fingerprints`
- `fingerprint_version`
- `fingerprint_algorithm`

The fingerprint representation is explicitly versioned so serialized objects remain interpretable if canonical serialization changes later.

Fingerprint version identifies the canonical serialization contract used to generate the digest. Fingerprint algorithm identifies the digest algorithm separately.

Fingerprints are internal integrity/audit aids, not substitutes for scientific provenance or source hashes.

Changing canonical fingerprint serialization requires a new `fingerprint_version`.

---

## 5. Feature representation

The molecule universe should avoid repeated per-molecule feature strings where practical.

Preferred representation:

- `features`: one-row-per-feature `DataFrame` containing `feature_id`
- `molecules$feature_index`: integer index into the feature census

This preserves feature identity while reducing repeated-string memory overhead.

---

## 6. `SegmentationPartition`

### Responsibility

One certified segmentation of one immutable molecule universe.

### Frozen conceptual slots

- `name`: `character`
- `cells`: `S4Vectors::DataFrame`
- `status`: `integer`
- `cell_index`: `integer`
- `ambiguity`: `S4Vectors::SimpleList`
- `run_metadata`: `S4Vectors::SimpleList`
- `provenance`: `S4Vectors::SimpleList`
- `universe_fingerprint`: `character`
- `assignment_fingerprint`: `character`
- `fingerprint_version`: `character`
- `fingerprint_algorithm`: `character`
- `schema_version`: `character`
- `metadata`: `S4Vectors::SimpleList`

### Cell census

`cells` has one row per segmentation cell and preserves original cell identifiers losslessly.

Required column:

- `cell_id`

The census includes zero-molecule cells and cells that later have zero correspondence.

### Assignment status

Frozen semantics:

- `assigned`
- `unassigned`
- `ambiguous`

Canonical v0.99 physical storage is an internal integer-coded vector.

For `assigned`, `cell_index` must reference a valid census row.

For `unassigned` and `ambiguous`, `cell_index` must be `NA_integer_`.

The integer codes are implementation details; public accessors expose semantic status names.

### Ambiguity

`ambiguity` contains:

- `evidence`: `DataFrame`
- `candidates`: `Hits`

`evidence` records ambiguous molecules and supports:

- candidate identities known
- candidate count known but identities unavailable
- candidate set entirely unknown

Suggested fields:

- `molecule_index`
- `candidate_count`
- `candidate_ids_known`

`candidates` stores only known molecule-to-candidate-cell relationships.

Ambiguous status with unavailable candidate identities is allowed and must never be fabricated or converted into unassigned status.

### Must not contain

- geometry
- expression matrices
- biological states
- correspondence edges
- support-filtered graphs
- pairwise comparison summaries

### Mutability

Conceptually immutable after construction. No public replacement methods for scientific core slots.

---

## 7. `SegmentationSet`

### Responsibility

One immutable molecule universe plus a named set of universe-bound segmentations.

### Frozen conceptual slots

- `molecules`: `S4Vectors::DataFrame`
- `features`: `S4Vectors::DataFrame`
- `segmentations`: `S4Vectors::SimpleList`
- `specimen_id`: `character`
- `coordinate_info`: `S4Vectors::SimpleList`
- `provenance`: `S4Vectors::SimpleList`
- `universe_fingerprint`: `character`
- `fingerprint_version`: `character`
- `fingerprint_algorithm`: `character`
- `schema_version`: `character`
- `metadata`: `S4Vectors::SimpleList`

### Molecules

Required scientific content includes:

- lossless `molecule_id`
- valid `feature_index`

Optional immutable columns may include:

- x/y/z coordinates
- filtering metadata
- other source molecule metadata

Coordinates are not required for correspondence-only workflows.

### Segmentations

`segmentations` is a uniquely named list of `SegmentationPartition` objects.

Every partition must match the set's molecule count, universe fingerprint, fingerprint version, and compatible fingerprint algorithm contract.

### Coordinate metadata

`coordinate_info` may declare:

- units
- origin
- axis directions
- transforms
- registration information

The slot records declarations; explicit scientific validation determines whether a spatial gate is satisfied.

### Must not contain

- `sf` geometry
- pairwise comparison results
- state-stability outputs
- dense correspondence matrices
- replicate-averaged results

---

## 8. `SegmentationComparison`

### Responsibility

Immutable complete pairwise correspondence between two segmentations before support thresholding, anchor selection, geometry restriction, or biological-state selection.

### Frozen conceptual slots

- `segmentation_names`: `character`
- `cell_summary_1`: `S4Vectors::DataFrame`
- `cell_summary_2`: `S4Vectors::DataFrame`
- `edges`: `S4Vectors::Hits`
- `components`: `S4Vectors::DataFrame`
- `assignment_summary`: `S4Vectors::DataFrame`
- `universe_size`: `numeric`
- `universe_fingerprint`: `character`
- `assignment_fingerprints`: `character`
- `fingerprint_version`: `character`
- `fingerprint_algorithm`: `character`
- `coordinate_info`: `S4Vectors::SimpleList`
- `provenance`: `S4Vectors::SimpleList`
- `schema_version`: `character`
- `metadata`: `S4Vectors::SimpleList`

### Cell summaries

Each side retains the full census and includes at least:

- `cell_id`
- `assigned_mass`
- `matched_mass`
- `unmatched_mass`
- `partner_count`
- `dominance_count`
- `component_id`

Zero-molecule and zero-correspondence cells remain represented.

### Canonical weighted graph

Canonical storage is `S4Vectors::Hits`.

- query indices refer to segmentation 1 census rows
- subject indices refer to segmentation 2 census rows
- only positive edges are stored

Edge metadata contains at least:

- `weight`
- `dominant_1`
- `dominant_2`
- `reciprocal_dominant`

Before object creation:

1. raw shared assignments are grouped by cell-index pair;
2. duplicate pairs are aggregated;
3. exactly one positive edge remains per pair;
4. canonical ordering is query index then subject index.

Original cell IDs are recovered through the side-specific census and are not copied onto every edge.

### Components

Component membership lives on cells through `component_id`.

The component table contains at least:

- `component_id`
- `n_segmentation_1`
- `n_segmentation_2`
- `component_class`

Classes:

- `one_to_one`
- `split`
- `merge`
- `many_to_many`
- `unmatched`

These must never be conflated with reciprocal-anchor `*_like` topology classes.

### Must not contain

- molecule-level assignment vectors
- full molecule universe
- geometry
- biological states
- thresholded replacement graphs
- anchor-only replacement graphs

---

## 9. Count and mass types

Identifiers and counts follow different rules.

Identifiers must never use ordinary double.

Scientific counts may use exact whole-valued R numeric/double values when they are:

- finite
- non-negative
- whole-valued
- `<= 2^53`

This applies to:

- edge weights
- assigned mass
- matched mass
- unmatched mass
- resolved mass
- same-state mass
- global molecule counts

If an exact count would exceed double's exact integer range, SegStability must fail rather than silently approximate.

---

## 10. Geometry ownership and optional sf dependency

Geometry is external to all three core classes.

`sf` belongs in `Suggests`.

No S4 slot requires an `sf` class, so assignment-only workflows can load and use SegStability without `sf` installed.

Spatial functions accept externally supplied geometry aligned by explicit segmentation cell IDs.

### Frozen spatial signatures

```r
anchorTopology(
    comparison,
    geometry,
    window,
    coordinate_info = NULL
)
```

```r
compareGeometry(
    comparison,
    geometry,
    coordinate_info = NULL
)
```

`geometry` is a named list keyed by the segmentation names in the comparison.

Geometry rows must expose explicit cell identity.

Unknown geometry cell IDs and duplicate geometry cell IDs are errors.

Missing geometry for census cells remains explicit geometry-unavailable status.

Invalid, empty, and non-positive-area geometries remain explicit and ineligible; they are not silently repaired.

### Frozen v0.99 observation-window representation

`anchorTopology()` uses a named numeric rectangular window:

```r
c(
    xmin = ...,
    xmax = ...,
    ymin = ...,
    ymax = ...
)
```

Requirements:

- numeric and finite
- exactly the names `xmin`, `xmax`, `ymin`, `ymax`
- `xmin < xmax`
- `ymin < ymax`
- expressed in the certified coordinate system used for geometry

Strict bbox interior eligibility uses this frozen rectangular observation window.

Arbitrary polygonal observation windows are deferred beyond v0.99.

---

## 11. Structural validity versus scientific validation

### Layer 1: `validObject()`

Validity checks inexpensive structural invariants only.

#### `SegmentationPartition`

- valid scalar name
- status/cell-index lengths agree
- only allowed status codes
- assigned entries have valid cell indices
- unassigned/ambiguous entries have `NA` cell index
- cell IDs structurally valid and unique
- ambiguity indices in range
- candidate indices in range
- ambiguity evidence compatible with status
- fingerprint fields structurally present

#### `SegmentationSet`

- molecule IDs structurally valid and unique
- feature census valid and unique
- feature indices in range
- segmentation names unique
- partition lengths match universe
- partition universe fingerprints match set
- fingerprint version/algorithm metadata structurally coherent

#### `SegmentationComparison`

- exactly two segmentation names
- Hits dimensions equal side census sizes
- positive whole edge weights
- unique edge pairs
- canonical edge order
- valid edge indices
- matched mass does not exceed assigned mass
- component IDs structurally valid
- fingerprint fields structurally coherent

### Layer 2: explicit scientific validation

Public validation functions:

- `validateSegmentation()`
- `validateSegmentationSet()`
- `validateCoordinates()`
- `validateGeometry()`
- `validateProvenance()`

These return inspectable validation reports and may support strict stop-on-failure behavior.

Scientific validation covers identity alignment, assignment conservation, coordinate compatibility, geometry eligibility, provenance completeness, and other gates that are more than cheap S4 structure checks.

No public `validation_level = "none"` escape hatch is provided in v0.99.

Constructors perform mandatory identity/alignment checks once. Analysis functions rely on constructed/fingerprinted objects plus structural validity and do not repeatedly redo expensive source-ID matching.

Spatial functions must run the coordinate/geometry gates necessary for safe spatial operation.

---

## 12. Provenance model

No dedicated provenance S4 class is introduced.

Provenance is a standardized `S4Vectors::SimpleList` with top-level concepts such as:

- `source`
- `source_version`
- `retrieval_date`
- `source_hashes`
- `software`
- `software_version`
- `software_hashes`
- `command`
- `parameters`
- `model`
- `threads`
- `runtime`
- `observation_window`
- `output_hashes`

Unknown values remain explicitly unknown, typically as typed `NA` or an empty structured value appropriate to the field.

Required provenance keys are not silently omitted solely because their values are unknown.

`validateProvenance()` derives completeness/status rather than relying on a mutable cached truth flag.

---

## 13. Biological-state API

States remain external to segmentation objects.

v0.99 supports categorical endpoints only.

Accepted user forms include:

- named character vector
- named factor
- named integer categorical labels
- `DataFrame(cell_id, state)`

State mapping must explicitly cover the requested full census.

Semantics:

- cell absent from mapping: malformed input/error
- cell present with `state = NA`: explicit unresolved state
- unknown cell ID in mapping: error
- duplicate state-map cell ID: error

Factors are interpreted by labels, not underlying integer codes.

Continuous endpoints are deferred.

Public function:

```r
stateStability(
    comparison,
    states
)
```

`states` is a named list keyed by compared segmentation names.

Return type: `S4Vectors::SimpleList` containing at least global, cell-level, and edge-level results.

---

## 14. Support sensitivity

Public function:

```r
supportSensitivity(
    comparison,
    thresholds = c(0, 0.01, 0.025, 0.05, 0.10, 0.20),
    direction = NULL,
    states = NULL
)
```

Return type: `S4Vectors::DataFrame`.

Thresholding never creates a new primary `SegmentationComparison` and never replaces:

- original denominators
- dominance
- components
- census

Thresholded edge selections are transient/internal views only.

Rows report support retained/lost together with conditional and total-population estimands as applicable.

---

## 15. Stochastic replicate analysis

Separate stochastic runs remain ordinary `SegmentationPartition` objects.

`run_metadata` may include:

- `method`
- `run_id`
- `replicate_group`
- `seed`
- `software_version`

Equal-valued cell IDs across runs have no matching semantics.

Public function:

```r
replicateAgreement(
    object,
    segmentations
)
```

where `object` is a `SegmentationSet` and `segmentations` explicitly names runs to compare.

Return type: `SimpleList` containing at least pairwise summaries and assignment-status contingencies.

Core results include:

- assignment-status agreement
- jointly uniquely assigned population
- label-invariant partition agreement such as ARI

Graph/state/topology variation remains explicit by analyzing each run separately rather than averaging results into one opaque statistic.

---

## 16. Frozen public API

### Constructors

- `SegmentationPartition()`
- `SegmentationSet()`

There is no public `SegmentationComparison()` constructor.

### Primary analysis

- `compareSegmentations()`

### Biological sensitivity

- `stateStability()`
- `supportSensitivity()`

### Spatial extension

- `anchorTopology()`
- `compareGeometry()`

### Replicate analysis

- `replicateAgreement()`

### Scientific validation

- `validateSegmentation()`
- `validateSegmentationSet()`
- `validateCoordinates()`
- `validateGeometry()`
- `validateProvenance()`

### Standard methods

- `show()`
- `summary()`

---

## 17. Frozen accessor API

- `segmentationNames()`
- `segmentationInfo()`
- `moleculeData()`
- `featureData()`
- `cellCensus()`
- `assignmentStatus()`
- `ambiguityData()`
- `correspondenceEdges()`
- `cellSummary()`
- `componentSummary()`
- `assignmentSummary()`
- `provenance()`
- `runMetadata()`

Avoid public generic names such as `edges()`, `edgeData()`, `components()`, `xSummary()`, or `aSummary()`.

Segmentation names, not mathematical X/A notation, are the preferred public selectors.

---

## 18. Comparison-engine rules

`compareSegmentations(object, segmentation1, segmentation2)` returns a complete `SegmentationComparison`.

The engine must:

1. operate on universe-aligned assignment vectors;
2. include only molecules uniquely assigned on both sides in positive-edge counting;
3. aggregate duplicate pair observations before edge creation;
4. compute each side's assigned denominator independently from its own assignment vector;
5. preserve zero-molecule/zero-correspondence census cells;
6. compute partner count from positive edges;
7. compute dominance from integer shared-molecule weight on the complete graph;
8. retain all exact maximum ties;
9. identify reciprocal dominance only for unique maxima on both endpoints;
10. compute connected components including isolates;
11. retain assignment-status cross-accounting;
12. never derive own-side denominators from edge row/column sums.

---

## 19. Memory/scaling strategy

v0.99 core objects are in-memory.

Approximate payloads excluding container overhead:

| Quantity | 1M | 10M | 50M |
|---|---:|---:|---:|
| one 32-bit integer vector | ~3.8 MiB | ~38 MiB | ~191 MiB |
| one 64-bit vector | ~7.6 MiB | ~76 MiB | ~381 MiB |
| status + cell-index per segmentation | ~7.6 MiB | ~76 MiB | ~381 MiB |
| x+y doubles | ~15 MiB | ~153 MiB | ~763 MiB |

Character IDs with all-unique values can dominate memory because each entry requires pointer storage plus string objects; integer/64-bit source IDs must therefore not be converted to unique strings merely for canonicalization.

Likely dominant terms at scale are:

1. source character molecule IDs when genuinely string-valued;
2. coordinates;
3. status/cell-index vectors across multiple segmentations;
4. retained molecule metadata;
5. geometry, which is deliberately external.

Correspondence storage scales with positive edges, not the Cartesian product of cell censuses.

Future backends may include HDF5/DelayedArray/Arrow-style storage, but none is required in v0.99.

The object model must not assume that assignment vectors can never become disk-backed later.

---

## 20. Interoperability scope for v0.99

### Required ingestion

- base `data.frame`
- `S4Vectors::DataFrame`

These require explicit scientific IDs; row order is never manufactured into molecule identity.

### Conditional interoperability

#### MoleculeExperiment

Supported only when original molecule identity is explicitly retained and can be proven to align with the SegStability universe.

If identity has been dropped by upstream processing, the adapter must fail.

#### SpatialFeatureExperiment

Useful as an optional geometry source when explicit cell identity and coordinate metadata are available.

Geometry presence alone does not establish molecule-assignment identity.

### Not promised as general SegmentationSet coercions

#### SpatialExperiment

Usually cell/spot-level and generally insufficient to reconstruct the immutable molecule universe.

#### SpatialExperimentIO outputs

Used according to the actual returned class and retained identities. SegStability does not duplicate platform importers or manufacture missing molecule identity.

Platform-specific raw readers are outside v0.99 scope.

---

## 21. Dependency plan for v0.99

### Imports

- `methods`
- `S4Vectors`
- `digest`

`S4Vectors` provides the core classes required by the object model. `digest` is used for versioned universe/assignment fingerprints.

### Suggests

Expected:

- `sf`
- `bit64`
- `MoleculeExperiment`
- `SpatialFeatureExperiment`
- `testthat`
- `BiocStyle`
- `knitr`
- `rmarkdown`

Packages should only remain in Suggests when executable code/tests/vignettes actually use them.

### Not dependencies initially

- `IRanges`
- `Matrix`
- `SparseArray`
- `DelayedArray`
- `SpatialExperimentIO`
- `igraph`
- `arrow`
- `HDF5Array`

No package is imported merely because it is common in Bioconductor.

---

## 22. Semantic immutability

R cannot prevent deliberate direct slot manipulation, so immutability is semantic rather than absolute.

Rules:

- no public replacement methods for scientific core slots;
- derived analyses return new tables/lists;
- support thresholding never modifies the primary graph;
- states are never written into `SegmentationComparison`;
- geometry is never written into core objects;
- fingerprints retain traceability to original aligned inputs;
- each class stores `schema_version` for future serialized-object migration.

---

## 23. Functions remaining internal

Examples:

- identifier-class validation
- LLint/integer64 exactness checks/helpers
- universe canonical serialization
- fingerprint generation
- external assignment reordering
- cell-ID-to-index translation
- grouped edge counting
- edge canonicalization
- denominator counting
- dominance/tie logic
- reciprocal-dominance logic
- component traversal
- state-map alignment
- support-edge selection
- contingency/ARI construction
- geometry spatial-index helpers
- geometry matching helpers

Low-level arithmetic functions are not exported merely for convenience.

---

## 24. Test obligations before substantive analysis code

The implementation must port the existing adversarial scientific oracle into package tests.

Particular preconditions before accepting LLint/integer64 support:

- exact duplicate detection
- exact matching
- exact ordering
- exact subsetting
- exact DataFrame round trip
- exact RDS serialization round trip
- no ordinary double conversion at any point

Constructor/alignment tests must cover reordered assignments, missing/extra IDs, duplicate IDs, unknown cells, and unsafe IDs before the comparison engine is considered trustworthy.

---

## 25. Deferred beyond v0.99

- continuous biological endpoints
- arbitrary polygonal observation windows
- disk-backed assignment classes
- Arrow/HDF5 assignment backends
- DelayedArray assignment backend
- frozen-reference training/clustering helpers
- platform-specific raw importers
- automatic segmentation-tool execution
- custom plotting framework
- thresholded primary-comparison objects
- consensus segmentation
- segmentation ranking/truth scores
- biological inferential hypothesis testing
- automatic geometry repair
- automatic coordinate transformation

---

## 26. Frozen workflow shape

### Geometry-free

```r
native <- SegmentationPartition(
    molecules = molecules,
    name = "native",
    cells = native_cells,
    molecule_id = native_assignment_ids,
    status = native_status,
    cell_id = native_cell_assignment,
    provenance = native_provenance
)

baysor <- SegmentationPartition(
    molecules = molecules,
    name = "baysor",
    cells = baysor_cells,
    molecule_id = baysor_assignment_ids,
    status = baysor_status,
    cell_id = baysor_cell_assignment,
    ambiguity = baysor_ambiguity,
    provenance = baysor_provenance
)

ss <- SegmentationSet(
    molecules = molecules,
    segmentations = list(native = native, baysor = baysor),
    specimen_id = "mouse_brain_roi1",
    coordinate_info = coordinate_info
)

cmp <- compareSegmentations(ss, "native", "baysor")

correspondenceEdges(cmp)
cellSummary(cmp, segmentation = "native")
componentSummary(cmp)
assignmentSummary(cmp)
```

### State sensitivity

```r
states <- list(
    native = DataFrame(cell_id = native_ids, state = native_state),
    baysor = DataFrame(cell_id = baysor_ids, state = baysor_state)
)

stab <- stateStability(cmp, states)
supp <- supportSensitivity(cmp, states = states)
```

### Optional geometry

```r
geometry <- list(native = native_sf, baysor = baysor_sf)
window <- c(xmin = 800, xmax = 1400, ymin = 500, ymax = 1100)

anchors <- anchorTopology(
    cmp,
    geometry = geometry,
    window = window,
    coordinate_info = geometry_coordinate_info
)

geom_cmp <- compareGeometry(
    cmp,
    geometry = geometry,
    coordinate_info = geometry_coordinate_info
)
```

### Stochastic replicates

```r
rep_agreement <- replicateAgreement(
    ss,
    c("proseg_1", "proseg_2", "proseg_3")
)
```

Run-local cell IDs are never matched literally across runs.

---

# FINAL FREEZE DECLARATION

**SEGSTABILITY_API_SPECIFICATION.md is FROZEN FOR v0.99 IMPLEMENTATION.**

The architecture may not be reopened merely for convenience during coding.

A change to this specification requires a concrete implementation contradiction, a failing scientific invariant that cannot be repaired within the frozen design, or a demonstrated Bioconductor compatibility requirement.

Ordinary implementation choices, optimization opportunities, naming preferences, or convenience requests do not constitute architecture contradictions.

The next phase is implementation planning and then package implementation against this specification.
