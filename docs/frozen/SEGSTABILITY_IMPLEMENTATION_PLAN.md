# SEGSTABILITY_IMPLEMENTATION_PLAN.md

Status: **FROZEN FOR v0.99 IMPLEMENTATION**

Companion document:
- `SEGSTABILITY_API_SPECIFICATION.md` — **FROZEN FOR v0.99 IMPLEMENTATION**

Purpose:
This document defines the implementation sequence for SegStability v0.99.x.
It does not redefine the scientific contract or public API. If implementation
reveals a concrete contradiction with the frozen API specification, that
specific contradiction must be documented and resolved before changing the API.

---

# 1. Implementation principles

1. Do not modify the frozen scientific-validation repository except for explicit
   documentation/export needs.
2. Create the package in a separate repository, expected at:
   `~/projects/SegStability`
3. Implement correctness and identity safety before optimization.
4. Port the existing adversarial scientific oracle into package tests as early
   as possible.
5. Do not implement geometry, state propagation, replicate analysis, adapters,
   or convenience helpers before the core universe/assignment/correspondence
   model is stable.
6. No dense cell-by-cell correspondence matrix is permitted as canonical storage.
7. No ordinary floating-point representation is permitted for scientific IDs.
8. No support-threshold operation may mutate or redefine the complete primary
   correspondence graph.
9. Geometry remains optional and must not become a package-load dependency.
10. No new scientific validation phase is required unless implementation exposes
    an actual contradiction in the frozen contract.

---

# 2. Target development environment

## Environment requirement

The implementation target is:

- the **current Bioconductor devel release at implementation time**;
- the **R version required by that Bioconductor devel release**;
- WSL/Linux development environment;
- Git repository separate from the prototype repository.

For reproducibility, the implementation record should note the versions expected
at the time this plan was frozen, but those versions are not a permanent
requirement. Before implementation begins, verify the current official
Bioconductor devel environment and use the R version required by that release.

At the time this plan was frozen (2026-10-06), the official Bioconductor install
page reports Bioconductor 3.24 as devel and R 4.6.0 as its required/supported R
version. This is a timestamped reproducibility record, not a permanent pin.

Because Bioconductor devel transitions with the R release cycle, revalidate the
development environment whenever the official devel requirements change and
rerun the package checks under the new devel environment.

## Environment gate

Before substantive implementation:

1. Verify the current official Bioconductor devel release and its required R version.
2. Install/use that R + Bioconductor devel environment.
3. Run:
   `BiocManager::valid()`
4. Confirm devel repositories are active and correct.
5. Record:
   - R version
   - Bioconductor version
   - configured repositories
   - platform
   - compiler
   - sessionInfo()

Package-name availability is a separate authoritative gate and must not be inferred
solely from whether `BiocManager::install("SegStability")` succeeds or fails.

## Authoritative package-name collision gate

Before package initialization:

1. Search the official current Bioconductor package namespace.
2. Search official historical/removed Bioconductor package records where available.
3. Search the official current CRAN package namespace.
4. Search official CRAN archive/historical package records.
5. Perform comparisons case-insensitively.
6. Record the search date, sources, and result.

The authoritative decision comes from these official current and historical
namespace searches.

`BiocManager::install("SegStability")`, if attempted, is only supplemental evidence.
An install failure must never be treated as proof that the package name is available,
because installation may fail for repository, network, dependency, version, or
environment reasons.

---

# 3. Repository/package initialization

Create:

`~/projects/SegStability`

Initial structure:

```text
SegStability/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── R/
├── man/
├── tests/
│   └── testthat/
├── vignettes/
├── inst/
│   └── extdata/        # only if tiny files are genuinely needed
├── NEWS.md
├── README.md
└── .github/
    └── workflows/
```

Initial package version:

`0.99.0`

Do not copy historical prototype code wholesale into the package.

Only port code after deciding whether each routine belongs in:
- production R code,
- tests,
- vignette support,
- or nowhere.

## Initialization gate

Package must:

- load successfully
- build successfully
- contain no substantive analytical implementation yet
- have valid DESCRIPTION/NAMESPACE scaffolding
- use version 0.99.0

---

# 4. DESCRIPTION strategy

## Initial Imports

Start with:

```text
methods
S4Vectors
digest
```

Each dependency must have a concrete purpose.

### methods
S4 class/method infrastructure.

### S4Vectors
Required for:
- DataFrame
- Hits
- SimpleList
- LLint

### digest
Used for versioned universe and assignment fingerprints.

## Initial Suggests

Start with:

```text
testthat
BiocStyle
knitr
rmarkdown
```

Add later only when implementation reaches the relevant module:

```text
sf
bit64
MoleculeExperiment
SpatialFeatureExperiment
```

## Do not add initially

Do not add merely speculatively:

```text
IRanges
Matrix
SparseArray
DelayedArray
igraph
arrow
HDF5Array
SpatialExperimentIO
```

A dependency should be introduced only when production code or executable tests
actually use it.

---

# 5. NAMESPACE strategy

Use roxygen2-generated NAMESPACE during development.

Principles:

1. Export only frozen public API.
2. Avoid wildcard imports.
3. Avoid generic names known to collide with graph/spatial ecosystems.
4. Reuse existing Bioconductor generics where appropriate.
5. Prefer package-specific accessors such as:
   - correspondenceEdges()
   - componentSummary()
   - assignmentSummary()
6. Do not export internal arithmetic helpers.
7. Do not define replacement methods for scientific core slots.

The initial NAMESPACE should remain very small.

---

# 6. Phase 1 — Lossless ID compatibility laboratory

This phase must occur before implementing scientific classes.

## Goal

Prove which identifier representations can be safely supported without
floating-point identity loss.

## Required ID classes

Test:

- integer
- character
- S4Vectors::LLint
- bit64::integer64, when bit64 is installed

## Required exactness tests

For LLint and integer64-like inputs, test:

1. duplicate detection
2. equality
3. matching
4. ordering
5. subsetting
6. DataFrame storage
7. DataFrame extraction
8. serialization with saveRDS/readRDS
9. round-trip preservation
10. values above 2^53
11. values near signed 64-bit boundaries where supported
12. no hidden coercion through ordinary double

## Rejection tests

Reject:

- ordinary numeric/double IDs
- logical
- factor
- complex
- raw
- missing IDs
- NaN/Inf
- unsupported arbitrary classes

## Continuous oracle porting in this phase

Immediately port the relevant historical portability/adversarial tests for:

- accepted/rejected ID types
- duplicate IDs
- missing IDs
- large exact integer IDs
- character/UUID-like IDs
- serialization safety

## Gate

LLint/integer64 support enters the official implementation only after all
required exactness tests pass.

No fallback through ordinary double is allowed.

---

# 7. Phase 2 — Core identity and fingerprint utilities

Implement internal utilities before S4 constructors depend on them.

## Required internal capabilities

- validate ID type
- validate missingness
- detect duplicate IDs exactly
- compare identity sets
- map external order to universe order
- detect missing IDs
- detect extra IDs
- map cell IDs to cell_index
- dictionary-code feature identity
- validate feature_index
- canonicalize assignment status input
- build ambiguity relationships
- calculate universe fingerprint
- calculate assignment fingerprint

## Fingerprint contract

Every fingerprint must retain:

```text
fingerprint_version
fingerprint_algorithm
fingerprint_value
```

Recommended initial conceptual version:

```text
fingerprint_version = "1"
```

The canonical serialization format used to generate the digest must itself be
versioned.

No serialized object may rely on an unversioned digest whose input encoding
could become ambiguous later.

## Fingerprint scope

Universe/assignment fingerprints are **internal integrity-binding mechanisms only**.

They are used to:

- bind partitions to the exact molecule universe;
- detect accidental mismatch or corruption;
- support serialization and reproducibility checks.

They do **not**:

- prove source authenticity;
- replace source hashes;
- replace provenance validation;
- establish scientific provenance.

Source provenance remains governed by the provenance contract and
`validateProvenance()`.

## Fingerprint tests

Test:

- deterministic repeatability
- molecule-order sensitivity
- molecule-ID sensitivity
- type sensitivity where scientifically relevant
- assignment sensitivity
- serialization round trip
- preserved fingerprint_version
- preserved fingerprint_algorithm

Continue porting historical identity/alignment oracle cases here, including
missing/extra IDs and reordered molecule assignments.

## Gate

No S4 constructor implementation proceeds until identity alignment and
fingerprint tests pass.

---

# 8. Phase 3 — SegmentationPartition class

Implement the first package-owned S4 class.

## Responsibilities

Represent:

- one certified segmentation
- one explicit cell census
- one assignment status per universe molecule
- one segmentation-local cell_index for assigned molecules
- ambiguity evidence
- run metadata
- provenance
- universe and assignment fingerprints

## Implementation order

1. class definition
2. cheap validity method
3. constructor
4. show()
5. accessors
6. unit tests

## Constructor behavior

The constructor must:

1. receive the molecule universe
2. validate assignment molecule IDs
3. detect duplicates
4. detect missing/extra molecule IDs
5. align exactly to universe identity
6. reorder once if needed
7. map assigned cell IDs to integer cell_index
8. retain ambiguity evidence
9. calculate fingerprints
10. drop duplicated assignment-side molecule ID storage after canonicalization

## Mandatory constructor tests

Port the corresponding historical adversarial/oracle cases immediately as part
of constructor development:

- correct ordered input
- reordered molecule assignments
- duplicate molecule ID
- missing molecule ID
- extra molecule ID
- unknown assigned cell
- duplicate census cell
- census integrity
- zero-cell census
- empty segmentation
- assigned with missing cell
- ambiguous with forced cell
- unassigned with forced cell
- assigned/unassigned/ambiguous semantics
- ambiguity with known candidate IDs
- ambiguity with only candidate count
- ambiguity with unknown candidate identities
- universe assignment conservation

## Gate

SegmentationPartition must pass all identity/alignment cases before
SegmentationSet implementation starts.

---

# 9. Phase 4 — SegmentationSet class

## Responsibilities

Represent:

- immutable molecule universe
- feature census
- named SegmentationPartition collection
- specimen identity
- coordinate declarations
- provenance
- universe fingerprint

## Implementation order

1. class definition
2. validity
3. constructor
4. accessors
5. show()
6. summary()
7. tests

## Constructor checks

- molecule IDs unique and valid
- feature IDs unique and valid
- feature indices valid
- segmentation names unique
- each partition has same universe size
- each partition has matching universe fingerprint
- specimen ID present
- coordinate_info structurally valid
- provenance structurally valid

## Continuous oracle porting in this phase

Port and run the historical cases covering:

- exact common-universe binding
- census integrity across segmentations
- missing/extra universe members
- reordered source assignments
- assignment conservation
- incompatible universe fingerprints

## Gate

A SegmentationSet containing multiple correctly aligned partitions must survive
serialization and reconstruct identical fingerprints before comparison work starts.

---

# 10. Phase 5 — SegmentationComparison class

Define the result class before implementing the computation engine.

## Responsibilities

Represent:

- two compared segmentations
- complete positive weighted correspondence graph
- independent denominators
- matched/unmatched mass
- partner counts
- dominance and ties
- reciprocal dominance
- connected components including isolates
- assignment-status summary
- source fingerprints

## Implementation order

1. class definition
2. validity method
3. internal-only constructor
4. accessors
5. show()
6. summary()

There is no exported public SegmentationComparison constructor.

## Gate

The class must be capable of representing:
- exact identity
- zero correspondence
- isolates
- split
- merge
- many-to-many
- dominance ties

before the production comparison engine is attached.

---

# 11. Phase 6 — Core correspondence engine

Implement `compareSegmentations()`.

## Computation order

1. retrieve canonical aligned assignment vectors
2. calculate N1 from all uniquely assigned molecules in segmentation 1
3. calculate N2 independently from segmentation 2
4. identify molecules assigned on both sides
5. group `(cell_index_1, cell_index_2)`
6. aggregate duplicate pair observations
7. store only positive edges
8. canonicalize edge ordering
9. create Hits
10. calculate matched mass per cell
11. calculate unmatched mass
12. calculate positive-edge degree
13. determine complete-graph integer-weight maxima
14. retain all exact ties
15. identify unique reciprocal dominance
16. compute connected components including isolates
17. construct assignment-status cross table
18. create SegmentationComparison

## Critical invariant

Never derive N1 or N2 from edge row/column sums.

They originate from each segmentation's complete uniquely assigned population.

## Performance rule

No `n_cells_1 × n_cells_2` dense allocation.

The computation should scale roughly with:

```text
|U| + number of positive edges + number of census cells
```

## Continuous oracle porting in this phase

Port the corresponding historical graph fixtures while implementing each graph
operation:

- exact identity
- zero correspondence
- pure split
- pure merge
- many-to-many
- isolates
- one-sided and both-sided dominance ties
- independent directional denominators
- matched mass
- unmatched mass
- connected components

These fixtures are development tests for the comparison engine, not a later
post-hoc validation step.

## Gate

Do not begin state or geometry modules until every exact graph fixture passes.

---

# 12. Phase 7 — Complete historical oracle port confirmation

By this point, relevant oracle cases have already been ported incrementally
alongside the ID, fingerprint, partition, set, and comparison implementations.

This phase is the explicit completeness audit confirming that the full historical
scientific/adversarial portability contract has been represented in package-native
tests with no material case omitted.

Recommended completeness review:

## Identity and census
- numeric/integer IDs
- large integer IDs
- character IDs
- UUID-like IDs
- unsafe ID rejection
- duplicates
- empty census
- unknown cells

## Assignment
- assigned
- unassigned
- ambiguous
- known ambiguity
- unknown ambiguity
- global conservation

## Graph
- identity
- zero correspondence
- split
- merge
- many-to-many
- one-sided dominance tie
- both-sided dominance tie
- isolates

## Counts
- independent denominators
- matched mass
- unmatched mass
- exact edge weights

## Components
- one-to-one
- split
- merge
- many-to-many
- unmatched

## Support-invariance precursor
Verify denominator/dominance/component values are stable and retrievable.

## First major package gate

Run:

```text
R CMD build
R CMD check
R CMD check --as-cran
BiocCheck
```

Do not treat optional geometry failures as blockers yet because the geometry
module does not exist at this stage.

---

# 13. Phase 8 — Scientific validation functions

Implement the explicit second validation layer:

```text
validateSegmentation()
validateSegmentationSet()
validateProvenance()
```

Later geometry phase adds:

```text
validateCoordinates()
validateGeometry()
```

## Validation report design

Return inspectable structured results rather than only TRUE/FALSE.

At minimum report:

- check
- status
- severity
- message

Expected statuses may include:

```text
PASS
FAIL
UNKNOWN
NOT_APPLICABLE
```

Strict mode may stop on failure.

Do not implement mutable validation caches in v0.99.

---

# 14. Phase 9 — Biological state module

Implement only categorical states in v0.99.

## State alignment rules

Supported input:

- named character vector
- named factor
- named integer categorical vector
- DataFrame(cell_id, state)

Rules:

```text
cell missing from mapping -> malformed input / error
cell explicitly state=NA -> unresolved
unknown cell ID -> error
duplicate cell ID -> error
```

## Implement stateStability()

Calculate:

- M
- R
- S
- S/R
- S/M
- S/sum(N1)
- S/sum(N2)
- M-R
- unmatched assigned mass each side
- ambiguous assignment mass
- per-cell Mi/Ri/Si
- Si/Ni
- Si/Ri
- mass-weighted totals
- equal-cell means over estimable cells

## Tests

Explicitly cover:

- equal states
- unequal states
- one missing endpoint
- both missing
- zero R
- zero M
- zero N
- zero-correspondence cells
- unresolved cells retained

## Gate

StateStability must reproduce hand-computable fixtures exactly before support sensitivity begins.

---

# 15. Phase 10 — Support sensitivity

Implement:

`supportSensitivity()`

## Requirements

Use the complete comparison as immutable source.

For each direction and threshold:

- w/N_original
- preserve original denominators
- preserve original dominance
- preserve original components
- preserve full census

Report:

- retained edges
- retained positive-degree cells
- retained mass
- retained fraction
- excluded matched mass
- unresolved mass
- conditional agreement
- total-population agreement

## Important design constraint

Do not return a new SegmentationComparison.

Do not expose a public API that can mistake the thresholded edge view for a new
primary graph.

## Gate

Existing adversarial test for unchanged denominators under thresholding must pass.

---

# 16. Phase 11 — Optional geometry module

Only at this stage add:

```text
sf
```

to Suggests.

## Implement first

`validateCoordinates()`

Then:

`validateGeometry()`

Then:

`anchorTopology()`

Then:

`compareGeometry()`

## Frozen window contract

v0.99 public window representation:

```r
c(
    xmin = ...,
    xmax = ...,
    ymin = ...,
    ymax = ...
)
```

Requirements:

- numeric
- finite
- exactly four named values
- xmin < xmax
- ymin < ymax
- interpreted in the already certified coordinate system

Arbitrary polygon observation windows are deferred.

## anchorTopology() implementation

Use:

- complete-graph unique reciprocal dominance
- valid positive-area geometry
- strict bbox interior
- both endpoints interior

Do not recompute dominance after support filtering.

## Geometry tests

- holes
- multipolygons
- invalid geometry
- empty geometry
- zero-area
- overlapping polygons
- boundary-inclusive PIP
- ROI-boundary touching
- missing geometry
- unknown geometry cell IDs
- duplicate geometry cell IDs
- wrong CRS
- units mismatch
- specimen mismatch
- incorrect geographic metadata on local coordinates

No implicit repair.

## Gate

The core package must still load and all non-spatial tests must still pass with sf absent.

---

# 17. Phase 12 — Replicate module

Implement:

`replicateAgreement()`

## Inputs

- SegmentationSet
- explicit character vector of segmentation run names

Use run_metadata:

- method
- run_id
- replicate_group
- seed
- software_version

## Required outputs

At minimum:

- pairwise run summary
- assignment-status contingency
- jointly uniquely assigned count/fraction
- label-invariant ARI

## Critical test

Permute all run-local cell IDs while leaving molecule assignments structurally identical.

ReplicateAgreement must return the same result.

This proves it is not matching cells by literal ID.

## Scope

Do not collapse topology/state variation into one number.

Higher-level run-specific graph/state/topology comparisons remain separate.

---

# 18. Phase 13 — Adapters

Adapters come after the core package works.

## Required v0.99 path

Support:

- data.frame
- S4Vectors::DataFrame

These remain the reference ingestion path.

## MoleculeExperiment

Add conditional adapter only when:

- original molecule identity is explicitly retained
- exact universe linkage is possible
- no identity is inferred from row order

If identity is absent:
fail.

## SpatialFeatureExperiment

Primary v0.99 role:

- optional geometry source
- possibly molecule source only when explicit original molecule IDs survive

Geometry does not prove assignment identity.

## SpatialExperiment

No general SegmentationSet coercion.

May provide cell-state annotations through explicit mapping.

## SpatialExperimentIO

Do not duplicate importers.

Use objects created by existing import infrastructure only when their retained
information satisfies SegStability's identity contract.

## Gate

Every adapter must have explicit failure tests for missing original molecule identity.

---

# 19. Phase 14 — Documentation and vignette

Create documentation after the API is exercised by production tests.

## Required documentation

- package-level documentation
- class documentation
- constructor documentation
- accessor documentation
- validation-function documentation
- analysis-function documentation
- geometry documentation
- replicate documentation

## Main vignette

Suggested title:

`SegStability: segmentation sensitivity from shared molecules`

Workflow:

1. define immutable molecule universe
2. create two tiny segmentation partitions
3. construct SegmentationSet
4. compare segmentations
5. inspect correspondence edges
6. inspect independent denominators
7. inspect zero-correspondence cells
8. inspect ties/components
9. supply categorical external states
10. run stateStability()
11. run supportSensitivity()
12. optional geometry
13. anchor topology
14. stochastic replicate agreement
15. interpretation boundaries
16. ecosystem positioning

## Synthetic vignette data

Keep tiny.

Target approximately:

- 20–50 molecules
- 5–8 cells per segmentation
- one ambiguous molecule
- one unmatched cell
- one dominance tie
- one non-one-to-one structure

Do not bundle validation ROIs.

---

# 20. Phase 15 — Performance profiling

After correctness is established:

Profile at increasing synthetic scales.

Suggested development-only scales:

```text
100k molecules
1M molecules
5M molecules
```

These are **resource-conditional targets**, not mandatory benchmark sizes.

Before each large profiling run:

1. inspect available physical RAM;
2. do not count swap as usable working memory;
3. estimate expected memory demand;
4. set a conservative memory ceiling;
5. monitor observed resource use;
6. stop or skip the benchmark if predicted or observed usage is unsafe.

The 5M-molecule profile is optional on constrained hardware.

Profiling should scale only as far as safely demonstrated. No implementation
milestone requires forcing a large run that exceeds the available physical-memory
envelope.

Measure:

- constructor peak memory
- universe alignment
- comparison construction
- edge aggregation
- component traversal
- stateStability
- supportSensitivity

Look specifically for:

- accidental character duplication
- repeated copies of assignment vectors
- dense intermediate matrices
- expensive DataFrame copies
- unnecessary sorting
- geometry accidentally loaded into core workflows

Do not introduce Arrow/HDF5 merely because memory is non-zero.

Only add a backend when profiling demonstrates a genuine need.

---

# 21. R CMD check / BiocCheck milestones

## Milestone A — Skeleton

Run:

```text
R CMD build
R CMD check
```

Expected:
clean package structure.

## Milestone B — Core classes

After SegmentationPartition/SegmentationSet:

```text
R CMD build
R CMD check
```

## Milestone C — Core comparison + adversarial oracle

Run:

```text
R CMD build
R CMD check
R CMD check --as-cran
BiocCheck
```

`R CMD check --as-cran` is advisory and useful for finding additional issues.
Its output must be interpreted rather than treated mechanically.

This is the first major implementation gate.

## Milestone D — State/support

Repeat all above.

## Milestone E — Geometry/replicate modules

Repeat all above:
- with sf installed
- with sf absent for geometry-free tests

## Milestone F — Adapters/vignette

Run:

```text
R CMD build
R CMD check
R CMD check --as-cran
BiocCheck
```

plus full vignette build.

## Milestone G — Submission readiness

Primary correctness is judged in the current Bioconductor devel environment using:

- `R CMD check`
- `BiocCheck`
- Bioconductor devel installation/testing

Quality interpretation:

- **ERROR** -> blocker
- **WARNING** -> blocker unless explicitly justified and accepted
- **NOTE** -> inspect and resolve when avoidable; environment-, network-, or
  platform-related notes may be documented rather than treated as intrinsic
  package failures

`R CMD check --as-cran` remains a useful advisory check, but external
`--as-cran` noise is not automatically an implementation failure.

Require:

- zero ERROR
- zero unjustified WARNING
- no avoidable intrinsic package NOTE
- BiocCheck clean or documented acceptable notes
- devel installation clean
- BiocManager::valid() clean
- examples executable
- vignettes executable
- no oversized bundled files
- no undocumented exports
- no unused dependencies
- package loads without sf
- package works from clean library
- CI passes

---

# 22. Recommended implementation order summary

The implementation sequence is frozen as:

```text
1. environment
2. repository/package skeleton
3. DESCRIPTION/NAMESPACE
4. ID exactness tests
5. universe/fingerprint utilities
6. SegmentationPartition
7. SegmentationSet
8. SegmentationComparison
9. compareSegmentations
10. portability/scientific oracle tests
11. explicit validation functions
12. stateStability
13. supportSensitivity
14. geometry validation
15. anchorTopology
16. compareGeometry
17. replicateAgreement
18. adapters
19. vignette/documentation
20. performance profiling
21. final R CMD check/BiocCheck hardening
```

---

# 23. Implementation gates

Do not advance if the current gate fails.

## Gate 1 — ID safety
Lossless identifier behavior proven.

## Gate 2 — Universe safety
Exact alignment, duplicate/missing/extra detection and versioned fingerprints proven.

## Gate 3 — Assignment safety
SegmentationPartition invariants proven.

## Gate 4 — Multi-segmentation safety
SegmentationSet universe binding proven.

## Gate 5 — Core graph correctness
Exact graph fixtures and independent denominators proven.

## Gate 6 — Scientific-contract portability
Existing adversarial cases ported and passing.

## Gate 7 — State/support correctness
All S/R/M and support accounting cases passing.

## Gate 8 — Spatial correctness
Geometry and rectangular-window rules passing with no silent repair.

## Gate 9 — Replicate correctness
Label-invariant replicate comparison proven.

## Gate 10 — Bioconductor readiness
Build/check/BiocCheck/vignette/CI clean.

---

# 24. What must not happen during implementation

Do not:

- modify the scientific definition to simplify coding
- derive denominators from edge sums
- silently coerce IDs through double
- invent missing molecule identity
- force ambiguous assignments
- silently repair geometry
- make sf mandatory
- create dense cell-by-cell matrices
- let support filtering replace the complete graph
- treat missing state labels as disagreements
- match stochastic-run cells by same-valued ID
- bundle large validation data
- add platform readers already provided elsewhere
- create additional public S4 classes without a concrete implementation contradiction
- reopen novelty or empirical validation merely because implementation begins

---

# 25. Definition of implementation-ready

SegStability is ready for substantive implementation when all of the following
are true:

- API specification is frozen
- this implementation plan is accepted
- authoritative current/historical CRAN/Bioconductor name-collision gate is recorded
- current official Bioconductor devel environment is valid
- package repository exists separately from prototype
- ID exactness test harness is ready
- no unresolved architecture contradiction remains

At that point, substantive implementation begins with the ID/universe layer,
not with comparison or geometry code.


# 26. Frozen implementation start boundary

This plan is **FROZEN FOR v0.99 IMPLEMENTATION**.

Do not create another implementation-planning phase.

Do not reopen API design unless implementation exposes a concrete contradiction
with the frozen specification.

Execution begins only in this order:

A. current Bioconductor devel environment gate  
B. authoritative current/historical CRAN/Bioconductor package-name collision check  
C. separate `SegStability` repository/package skeleton  
D. lossless-ID exactness laboratory  
E. universe/fingerprint utilities

Do not begin:

- `compareSegmentations()`
- `stateStability()`
- geometry modules
- replicate analysis
- adapters
- vignettes

until their prerequisite gates in this plan have been reached.
