# SegStability final pre-API evidence audit

Date: 2026-10-05. Scope: the accepted mouse-brain source and two prospectively
selected, spatially disjoint ROIs; independent re-audit, replication and final
scientific contract. This supersedes the provisional Phase 7 API prerequisites,
not its immutable observations. Historical breast geometry-conditioned results
are not validation evidence. No package API was designed or implemented.

## 1. Executive summary

**Final recommendation: GO TO API DESIGN**, for a narrowly scoped descriptive
segmentation-sensitivity analytical layer, not a segmentation-quality truth score.
Both clean ROIs support the correspondence method and descriptive state-sensitivity
estimands. They do not establish biological generalization, universal topology
penalties, predictive superiority over geometry or new bipartite-graph mathematics.

ROI 2 completed Baysor and three separately preserved ProSeg executions on this
host. Its strict all-248-reference-feature attempt failed: Slc6a3 is absent in
eligible native cells. With explicit user approval, the already fixed ≥3-native-cell
feature rule now selects the fitted reference features. ROI 2 fits 247 genes;
normalization libraries, input feature order, assignments and graph denominators
still retain all 248 input features and all 301,784 evaluation molecules. The failed
attempt remains preserved. No segmentation or scientific parameters were tuned.

Independent reconstruction passed 49,741 checks for ROI 1 and 51,604 for ROI 2,
plus 229 additional checks of replicate ARI, geometric calculations and the exact
common-reference sensitivity. These are coverage evidence, not proof by check count.
Independent sf samples agree on every assignment tested after lossless coordinate
transport. Thirty-nine general scientific-contract adversarial cases passed.

Most qualitative conclusions replicate. Baysor native-direction split-like and
complex anchor strata have lower conditional stability in both ROIs; native-direction
merge-like strata have higher stability. Rare ProSeg contrasts vary in sign across
runs, ROI and reference. Anchor optimism is possible, not universal. Novelty is
**CONDITIONAL**: no complete integrated A–H implementation was identified in the
bounded maintained-source review, but clustComp already implements shared-element
weighted bipartite comparison. The defensible contribution is integrated scientific
accounting and propagation, not inventing that representation.

Machine-readable evidence is in [results/final_pre_api](results/final_pre_api/).
The final contract is [SEGSTABILITY_SCIENTIFIC_SPECIFICATION.md](SEGSTABILITY_SCIENTIFIC_SPECIFICATION.md).
The closure ledger and verification are `artifact_manifest.json`,
`artifact_hashes.sha256`, `final_gate_validation.json` and `final_run_manifest.json`
in that evidence directory. Final classifications are in section 19.

## 2. Independent Phase 7 re-audit

`141_independently_reaudit_clean_evidence.py` uses independent integer counters,
maxima, BFS components, PCA algebra and brute-force distances, not production
graph/stability/voting routines. `142` exports frozen R components and repeats
fixed-seed Leiden. `146/147` independently reconstruct sf geometries and sample
molecules deterministically by original transcript identity. `155` independently
calculates contingency ARI, scalar polygon areas and common-reference votes.

| Evidence tested | Independent result / limitation |
|---|---|
| Prospective ROI freeze | Frozen specification/hash precedes alternative generation; no state/topology-based selection |
| Sources/software | Accepted source hashes, archives and executable hashes verified; historical Phase 7 ledger rechecked at closure |
| Baysor | Frozen 0.7.1 command, no native-prior argument, same scientific settings |
| ProSeg | Three distinct run receipts/output hashes; native-prior dependence retained and declared |
| Assignment | Original evaluation transcript/native IDs unchanged; inclusive unique/zero/multiple-hit policy preserved; mass conserved |
| Edges/denominators | Rebuilt from transcript identity; own-side assigned denominators independently counted, not shared row sums |
| Dominance/topology | Complete-graph maxima, exact ties, unique reciprocity and interior fan-in/out independently checked |
| Components | Independent BFS includes isolates; numbers arbitrary, memberships/structure consistent |
| Native reference | Native-only eligibility, means/sample SD, PCA eigen-equations/projection, exact kNN/sNN and repeated Leiden consistent |
| Transfer | Independent exact distances and voting for every method, three resolutions and three k settings; unresolved cases retained |
| Primary/anchor/support | Recomputed S/R/M and own-side totals, equal-cell means, selection gaps and all six support thresholds |
| Replicates/geometry | Independent ARI arithmetic, IoU/Dice/areas and maximum-IoU objective agree; GEOS dependency is shared |

Numerical comparison uses absolute/relative tolerance 1e-10 for independently
recomputed floating summaries; integer identity/mass assertions are exact. Frozen
Leiden memberships repeat exactly in the recorded environment, not across arbitrary
future package versions. Geometry cross-check: 8,112 ROI 1 and 8,193 ROI 2 molecules,
each against all four alternatives. This is a sampled independent implementation
check alongside full production conservation, not an independent full reassignment.

Falsification found a genuine technical portability defect in the unchanged
reference recipe, a misleading generic anchor column name, reference-dependent
rare effects and limits on all-population topology interpretation. None invalidates
the independently recomputed full-graph mass estimands after the documented fix.
No CRITICAL clean-evidence defect was confirmed. Section 16 preserves failures and
qualifications instead of interpreting passing checks as universal certification.

## 3. ROI 2 selection and resource rationale

Rule registered before counting proposed windows: translate the ROI 1 midpoint
exactly 1,000 µm left, then select the smallest integer square half-width reaching
300,000 filtered molecules. Only XY, dimensions, molecule/native-cell counts,
boundary bounds and the measured Phase 7 resource envelope enter selection.
No topology, reference composition, agreement or biological output is inspected.

| Quantity | ROI 2 |
|---|---:|
| x interval (µm) | 1438.25–2010.25 |
| y interval (µm) | 1483.275–2055.275 |
| Dimensions / area | 572 × 572 µm / 327,184 µm² |
| Evaluation molecules | 301,784 |
| Native census / molecule-bearing cells | 1,077 / 1,072 |
| Native strict bbox-interior / crop-edge | 925 / 152 |
| Gap to ROI 1 in x | 425 µm |
| Specification frozen (UTC) | 2026-10-05 20:43:09.842468 |
| Baysor began (UTC) | 2026-10-05 20:43:24.466749 |

Actual host physical RAM was 13,603,368,960 bytes (12.67 GiB), available physical
RAM 8,703,246,336 bytes (8.11 GiB), free swap 16 GiB. Swap was not credited as
capacity. The budget was min(8 GiB, 70% available, available−3 GiB), giving
5,482,020,864 bytes (5.11 GiB). Output is on ext4; WSL2 Linux, i7-9750H,
6 cores/12 logical CPUs, approximately 796 GiB free disk; effective cgroup/ulimit
values are saved in `roi2/host_resource_manifest.json`.

The existing conservative envelope predicted 5,450,588,299 bytes and 206.74 seconds.
That prediction is just below the process ceiling, not a large ceiling margin;
the live guard was therefore essential. Actual peak was 4.758 GiB, runtime 210.58 s,
with no process swap. No repeated full calibration or failed full-subset job was run.
Monitoring retained resource checkpoints, commands and termination guards.

## 4. ROI 2 source/integrity gate

Dataset: `Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP`, the Phase 5 accepted
official source, XOA 1.0.2, research/example reuse under the recorded CC BY 4.0
license. Phase 5 source URLs, archive/version/CRC metadata and twelve input SHA-256s
remain authoritative and were reverified, not replaced. See
`results/fix_phase5/source_artifact_manifest.csv` and the ROI gate receipt.

Filtering is QV≥20 and Gene Expression features. Original transcript IDs, cell IDs
and native coordinates are preserved. Native unassigned is −1; ProSeg cell 0 is
a valid segmentation-specific ID, not an unassigned sentinel. Native assigned
molecules are 294,506; unassigned 7,278. Their sum is 301,784. Five census cells
receive no evaluation molecules. Cells, official feature matrix and boundary
namespace/coordinates pass the accepted-source/native-ID consistency gates.
Full-cell official counts include outside-window molecules and are not required
to equal ROI counts. Native assignments use exact mask IDs, not visualization PIP.

| Frozen object | SHA-256 |
|---|---|
| ROI 2 specification | `daae64d5940a41326bb1793ea5f4896b58631ddf0d695219062b189741db329d` |
| ROI 2 molecule Parquet | `c189fc2537c2db28fd763f03b80fb6ce769a33ba7b4115ea9109b236872a546f` |
| Row-equivalent Baysor CSV | `e00a48e8f7518c28281f34e7a44f8b35c289e49a4107207025688926d4851a7a` |
| ROI 1 frozen reference | `87a2c3f7a194c7c9902ea31ea6a78a6f55ec614af45de48b7bf51b66baf8b0e3` |
| Approved ROI 2 frozen reference | `c466fd62675ddca163d9b6c10ca3d0c4f437864469288cf5485b5acf123a2f51` |

CSV transport was read back under the original Arrow schema and compared exactly.
Later independent sf diagnostics use float32 coordinates promoted losslessly to
float64 and serialized with 17 significant digits; see section 16.

## 5. ROI 2 Baysor/ProSeg execution

Scientific Baysor command, changing only paths/window:

```text
baysor run -x x_location -y y_location -z z_location -g feature_name
  -m 20 --n-clusters=8 -s 10 --scale-std=25%
  --polygon-format=FeatureCollection -o OUTPUT molecules.csv
```

Baysor 0.7.1, eight Julia threads; no segmentation prior. Source/archive and binary
hashes are in the final run manifest; executable SHA-256 is
`db29893013af7b34a5d962132348a25d51b4bec87f342ea0c06d9bb86f97f791`.
Julia 1.10.6 bundle. Successful sampling alone did not qualify the run: assignments,
cell statistics, polygons and complete logs exported and were validated.

ProSeg 3.2.0 uses the identical native-prior design and separate stochastic runs:

```text
proseg --xenium RAW_SCHEMA.parquet --output-spatialdata spatialdata.zarr
  --output-cell-metadata cell_metadata.csv.gz
  --output-transcript-metadata transcript_metadata.csv.gz
  --output-cell-polygons cell_polygons.geojson.gz --nthreads 8
```

Executable SHA-256:
`846681b54975343b48f0cac829f6d8445e767f4cd4e87517e32a16e095ea2483`.
Three independent executions, no exposed fixed seed, are retained without averaging.
The raw-schema adapter has 379,522 spatial rows and preserves official columns/IDs;
the common scientific evaluation uses only the 301,784 qualified molecules.
Raw control-feature vocabulary differs from the 248 Gene Expression evaluation
features; this is intentional, unchanged from Phase 7 and recorded.

| Execution | Seconds | Peak RSS GiB | Geometry cells | Export gate |
|---|---:|---:|---:|---|
| Baysor | 210.58 | 4.758 | 1,203 | PASS |
| ProSeg 1 | 138.47 | approximately 0.207 | 1,016 | PASS |
| ProSeg 2 | 133.44 | approximately 0.207 | 1,016 | PASS |
| ProSeg 3 | 135.48 | approximately 0.207 | 1,016 | PASS |

No sustained process swap was observed. Full commands, environments, timestamps,
checkpoints, log/output hashes are in `roi2/alternative_segmentation_runs.json`.
ProSeg observed XY agrees with native coordinates; its transformed observed Z is
source-defined tilt/quantile-normalized, not micrometres. Certification tests that
transformation explicitly. Common PIP uses original native XY, never relocated
latent positions. No geometry buffer, smoothing, repair or undocumented import.

## 6. ROI 2 correspondence and structural results

All 1,077 native visualization polygons and all alternative polygons are valid,
nonempty and positive-area. Native/Baysor each have one part per cell; ProSeg runs
contain 1,286/1,282/1,308 parts in 1,016 objects, with multipart objects retained.
No holes are observed, though hole behavior is tested adversarially. Positive-area
overlap pairs are native 475, Baysor 1,305 and ProSeg 94/111/101. This is why ambiguity
must remain explicit. Native assignments are authoritative mask assignments; native
visualization overlaps do not replace them. Baysor polygons all lie within the ROI,
which does not prove freedom from window effects.

| Method | Assigned | Outside | Ambiguous | Shared mass M | Edges | Unique reciprocal edges | Interior anchors | Zero-edge X / A |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Baysor | 282,615 | 11,091 | 8,078 | 275,977 | 3,633 | 676 | 592 | 17 / 37 |
| ProSeg 1 | 210,607 | 90,857 | 320 | 210,542 | 2,238 | 1,006 | 915 | 51 / 1 |
| ProSeg 2 | 209,867 | 91,595 | 322 | 209,811 | 2,233 | 1,006 | 914 | 50 / 0 |
| ProSeg 3 | 209,174 | 92,316 | 294 | 209,112 | 2,230 | 1,006 | 915 | 50 / 0 |

Each assignment row sums to 301,784. Full census nodes, zero edges and ties remain
in `graphs/`; shared counts derive from common transcript identity and independent
own-side denominators. Complete-graph dominance is calculated before support or
anchor selection. Connected components include isolates. Weak shared edges can
produce a giant many-to-many component: a component label is not an unselected
equivalent of the four anchor-topology classes. Population flow, distributions,
geometry status, state composition and zero/undefined strata are in
`roi2_qualified/population_geometry/` and `stability/`.

## 7. ROI 2 biological sensitivity

The failed strict reference attempt is `roi2/replication_execution.json`. It stopped
at the all-248-feature ≥3-native-cell assertion, before any primary endpoint existed.
The approved amendment is `roi2_qualified/feature_policy_amendment.json`, frozen
before amended fitting. Slc6a3 has zero expression in eligible native cells; it is
excluded only from fitted PCA features. All 248 features still define library
sizes and all molecules still define assignments/edges/denominators.

One native-only reference uses 1,057 eligible native cells, 247 fitted features,
log2(1+5000×count/full-panel-library), native means/sample SD, full PCA rank 30,
PCs 1–15, exact k=30 sNN and modularity Leiden resolution .15/100 iterations/seed
2026. Resolution .10/.15/.30 produces 3/4/7 states (ROI 1: 4/7/9). All alternatives
use exactly that same frozen reference within ROI 2. One-way transfer uses k=30,
with k=15/50 sensitivities; ties, cutoff-distance ties and low-feature queries remain
unresolved. No posterior confidence threshold was fitted.

Primary resolved transfer confidence q10/median/q90 for Baysor is .6333/.9667/1;
ProSeg 1 q10 .7667 and ProSeg 2 .7333, medians 1. Confidence is vote concentration,
not a calibrated probability. Baysor has 65 low-feature cells (123 molecules) and
7 tied cells (1,802 molecules); ProSeg 1 has 1 low-feature cell (0 molecules), 4 tied
(282 molecules); ProSeg 2 has 1 low-feature (3), 3 tied (583). All run/settings/status
counts, confidence/margin/distance distributions are retained in
`transfer/transfer_diagnostics.csv`; endpoint-unresolved matched mass also includes
native unavailable states and is therefore not the same as query unresolved mass.

Let S be same-state shared mass, R resolved-label shared mass, M all shared mass;
N_X and N_A are independent uniquely assigned totals. These are primary complete-
population estimands, not anchor agreement:

| ROI / method | S | M−R | S/R | S/N_X | S/N_A |
|---|---:|---:|---:|---:|---:|
| 1 Baysor | 222,274 | 1,475 | .806741 | .748977 | .789954 |
| 1 ProSeg 1 | 189,895 | 214 | .927240 | .639873 | .925992 |
| 1 ProSeg 2 | 189,900 | 761 | .931535 | .639889 | .927817 |
| 1 ProSeg 3 | 188,960 | 1,556 | .929057 | .636722 | .921842 |
| 2 Baysor | 234,159 | 1,936 | .854467 | .795091 | .828544 |
| 2 ProSeg 1 | 197,936 | 285 | .941400 | .672095 | .939836 |
| 2 ProSeg 2 | 196,280 | 589 | .938142 | .666472 | .935259 |
| 2 ProSeg 3 | 195,449 | 228 | .935682 | .663650 | .934385 |

For example ROI 2 Baysor R=274,041, M=275,977; ProSeg 1 R=210,257, M=210,542.
S/M, both unmatched totals, ambiguous mass, equal-cell summaries and non-estimable
cells are explicit in `cross_roi/endpoints.csv` and the underlying per-cell tables.
Baysor exceeds ProSeg on native-assigned S/N_X but is lower on S/R and alternative-
assigned S/N_A in both ROIs. This is an estimand-dependent description, not an
algorithm-accuracy ranking or a causal native-prior effect.

## 8. Cross-ROI replication

The ROIs are disjoint parts of one specimen, not biological replicates. Local
reference-state IDs/granularity are not aligned labels. Values are not averaged.

| Qualitative claim | Cross-ROI adjudication |
|---|---|
| Alternatives create non-equivalent partitions | REPLICATED |
| Shared-molecule graph captures many-to-many relationships | REPLICATED |
| Full-population and anchor estimands differ | REPLICATED |
| Anchor optimism possible but not universal | PARTIALLY REPLICATED |
| Non-one-to-one structure does not uniformly lower state stability | REPLICATED |
| Some topology associations occur but depend on method | REPLICATED, selected/descriptive scope |
| Directional estimands matter | REPLICATED |
| Rare ProSeg topology effects are less stable than global endpoints | REPLICATED |
| Molecule graph retains information distinct from geometric matching | REPLICATED, not predictive-superiority evidence |

Important exception: all ROI 2 local-reference primary anchor gaps are positive,
whereas ROI 1 ProSeg run 3 has negative gaps. It would be incorrect to claim every
ROI demonstrates anchor pessimism. The common-reference sensitivity recovers a
small negative ROI 2 alternative-side run 3 gap, demonstrating reference dependence,
not replacing the primary result.

Secondary sensitivity projects BOTH ROI 2 native cells and alternatives into the
exact unchanged 248-gene ROI 1 reference. No retraining or state-ID alignment.
`roi2_common_reference/` is isolated and labelled secondary. Its S/R is Baysor
.881515 and ProSeg .956268/.954196/.957073. Native Baysor split/complex lower and
merge higher patterns survive; some tiny alternative-side merge and rare ProSeg
contrasts change sign. This separates some reference-granularity effects from
spatial robustness. Additional independent brute-force checks verify all native
queries and alternative votes for all nine settings in this branch.

## 9. Anchor-selection comparison

Three quantities must not be confused: unweighted reference-label agreement on
resolved anchor PAIRS, molecule-weighted stability of all edges incident to anchor
CELLS on each side, and molecule-weighted agreement restricted to anchor-pair EDGES.
`cross_roi/explicit_anchor_pair_agreement.csv` separates pair quantities; `anchors.csv`
contains side-specific incident-cell selection estimands. Historical generic
`resolved_conditional_stability` naming is qualified, not overwritten.

| ROI / method | Native anchor−nonanchor S/R | Alternative anchor−nonanchor S/R |
|---|---:|---:|
| 1 Baysor | +.082447 | +.246538 |
| 1 ProSeg 1 | +.017702 | +.001359 |
| 1 ProSeg 2 | +.017301 | +.011916 |
| 1 ProSeg 3 | −.001364 | −.010314 |
| 2 Baysor | +.029776 | +.223345 |
| 2 ProSeg 1 | +.029693 | +.023631 |
| 2 ProSeg 2 | +.029598 | +.020557 |
| 2 ProSeg 3 | +.010991 | +.001623 |

ROI 1 Baysor anchors represent 61.1% of native cells and 76.81% of native-incident
matched mass; ROI 2 54.97% and 69.22%. ROI 2 alternative-side anchors represent
49.21% of Baysor cells and 80.02% of incident matched mass. These are selected
populations. Example ROI 1 Baysor pair agreement .914614 versus weighted pair-edge
agreement .944123 is not the complete population .806741 or native anchor-cell
conditional stability. None may substitute for another.

## 10. ProSeg stochastic robustness

Three executions per ROI are shown separately. Same-valued run-specific cell IDs
are not matched literally; transcript identity defines replicate contingencies.
Cell counts are constant within ROI (944 in ROI 1, 1,016 in ROI 2), but partitions
are not identical. Jointly uniquely assigned proportions are .6483–.6490 in ROI 1
and .6650–.6676 in ROI 2. Conditional ARI is .960626–.962140 and .963787–.964798,
respectively. Status agreement over ALL molecules is .9372–.9397 and .9414–.9419.
High conditional ARI does not certify the outside/ambiguous population.

ROI 2 primary S/R spans .005718 and native S/N_X spans .008445. Rare native merge
strata (7/7/6 anchors) contrast with one-to-one by +.044577/−.023518/+.056328,
versus ROI 1 (6/6/5) −.475409/−.266776/+.067271 approximately. Those rare effects
are not a repeatable general penalty. ROI 2 split counts 1/3/3 and one complex
anchor only in run 3 further limit such claims. Full graph/assignment, contrast
and endpoint variation remain separate in `replicate_robustness/` and cross-ROI CSVs.

## 11. Geometry-vs-molecule comparison

Maximum-total-IoU positive one-to-one matching of visualization polygons omits
2,994 of 3,853 Baysor molecule edges in ROI 1 and 2,780 of 3,633 in ROI 2. ProSeg
one-to-one geometry matching omits 1,084–1,104 and 1,215–1,223 molecule edges,
respectively. Some positive geometry pairs have no shared molecules; some shared
edges have zero positive-area intersection. Boundary-inclusive points, visualization
approximation and exact native masks explain why equality is not required.

Scalar area/IoU/Dice recalculation agrees. Independent minimization of negative
IoU reproduces the saved maximum objective; it uses the saved all-positive-pair
table, not independent enumeration of every non-edge. Assignment samples use sf;
both geometric paths still share GEOS. The defensible result is additional
many-to-many/mass information beyond one-to-one matching, not new geometry
mathematics, a validated accuracy score or improved prediction.

## 12. Consolidated scientific specification

[SEGSTABILITY_SCIENTIFIC_SPECIFICATION.md](SEGSTABILITY_SCIENTIFIC_SPECIFICATION.md)
defines the lossless molecule universe/cell censuses, coordinate/provenance gates,
authoritative versus polygon assignment, status policy, independent denominators,
complete graph, dominance/ties/components, exact selected topology, frozen native
reference and transfer, all-edge/all-cell versus secondary anchor estimands, support
loss and inference/dependence boundaries. It includes the approved native-only
feature policy, full-panel normalization and exact unresolved-state handling.
This is a scientific contract, not an API or proposed S4 class.

All six thresholds (0,1,2.5,5,10,20%) and nine reference/transfer settings are saved
for full-observed, interior and anchor-pair support in both directions. At 20%, ROI 2
Baysor native-direction support retains 86.85% matched mass, excludes 36,282
molecules and raises conditional stability to .908982 while reducing S/N_X from
.795091 to .734980. Alternative direction retains 85.63%, excludes 39,652 and has
.892691 conditional / .741917 total. ProSeg native support retains approximately
95.2–95.4%, with S/N_X reduced to .642737–.651623. Conditional increases are never
reported without original-denominator support costs. Exact tables: `cross_roi/support.csv`.

## 13. Novelty landscape

Methods: executing-plans and literature-review skills guided implementation and
source logging. This is a one-reviewer focused rapid/scoping review dated 2026-10-05,
not an exhaustive systematic review or installed competitor benchmark. Queries,
scope, pinned source manifests and feature matrix are machine-readable in
`novelty_adjudication/`. Ten repository snapshots, full Giotto/GiottoClass sources
and official clustComp/clusterExperiment sources provide 395 file observations
(some Giotto files occur in both snapshots), not 395 independent methods.

A=multiple partitions; B=shared-identity weighted bipartite graph; C=independent
directional assignment denominators; D=many-to-many/dominance/topology; E=full-
population state sensitivity; F=unresolved/unmatched/support; G=stochastic partition
diagnostics; H=anchor-selection accounting. Y=explicit capability (sometimes generic),
P=partial/adjacent/derivable overlap, I=infrastructure, ND=integrated capability not
identified in inspected sources, U=unverified. ND does not prove absence.

| Candidate | A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|---|
| Sopa | Y | P | ND | P | P | P | P | ND |
| ResolVI/scVIVA-tools | P | ND | ND | P | P | P | P | ND |
| cellAdmix | P | P | P | P | P | P | P | ND |
| CRISP | P | ND | P | ND | P | P | P | ND |
| MARC | Y | U | U | P | P | P | U | U |
| MoleculeExperiment | Y | I | I | I | I | P | I | ND |
| SpatialExperiment | I | I | I | I | I | I | I | I |
| SpatialFeatureExperiment | I | I | I | I | I | I | I | I |
| SpatialExperimentIO | I | I | I | I | I | I | I | I |
| MultiAssaySpatialExperiment | Y | I | I | I | I | P | I | ND |
| Giotto/GiottoClass | Y | P | P | P | P | P | P | ND |
| clustComp | Y | Y | P | P | ND | P | P | ND |
| clusterExperiment | Y | P | P | P | P | P | Y | ND |

Sopa supports segmentation-layer comparison and transcript aggregation; its
geometric overlay is reconciliation, not this full mass/anchor contract.
[Official comparison tutorial](https://prism-oncology.github.io/sopa/tutorials/compare_segmentations/).
ResolVI/scVIVA provides generative denoising/state transfer; cellAdmix estimates
contamination/ownership, and CRISP evaluates marker coexpression. Their partial
overlap is substantive, not dismissed by differing names.
[scVIVA source](https://github.com/YosefLab/scviva-tools),
[cellAdmix source](https://github.com/kharchenkolab/cellAdmix-core),
[CRISP source](https://github.com/santangelo-lab/CRISP).
MARC is current consensus-segmentation prior art; a public maintained implementation
was not located, so its code capabilities remain unverified rather than absent.
[MARC preprint](https://arxiv.org/abs/2609.13665).

Most importantly, clustComp explicitly compares clusterings using intersection-
weighted bipartite graphs; inspected `flatVSflat` builds `table(flat1,flat2)`.
The graph core is established prior art. clusterExperiment adds multiple-clustering,
unassigned-label and repeat-clustering infrastructure.
[clustComp](https://bioconductor.org/packages/release/bioc/html/clustComp.html),
[clusterExperiment](https://bioconductor.org/packages/release/bioc/html/clusterExperiment.html).

Verdict: **NOVELTY CONDITIONAL** for the integrated A–H spatial sensitivity workflow.
No complete implementation was identified in bounded sources; no claim of being
first at segmentation stability, molecule counting, partition graphs or label
transfer. Scope/support accounting and full-population propagation are the potential
software gap. This qualified contribution can be specified for API design without
another generic novelty phase; publication priority is not proven.

## 14. Bioconductor fit

Recommend **NEW PACKAGE** as a small analytical layer reusing established experiment
and geometry infrastructure, not a new universal spatial-data container or wrapper
that requires Baysor/ProSeg execution to compare externally certified assignments.
The combination of accounting, state propagation and selection diagnostics is
larger than a demonstration vignette, while the input models already exist.

MoleculeExperiment represents molecules/boundaries and counts into SpatialExperiment;
its essential-column import can drop transcript identity, so identity must be
explicitly retained. Overlapping-polygon counting does not implement the required
unique/ambiguous assignment contract automatically.
[MoleculeExperiment manual](https://bioconductor.org/packages/release/bioc/manuals/MoleculeExperiment/man/MoleculeExperiment.pdf).
SpatialExperiment/SummarizedExperiment cell columns are not interchangeable
segmentations with different cell universes; SpatialFeatureExperiment supplies
geometry/spatial infrastructure, and SpatialExperimentIO import utilities.
[OSTA infrastructure](https://bioconductor.org/books/release/OSTA/pages/bkg-infrastructure.html),
[SpatialFeatureExperiment source](https://github.com/pachterlab/SpatialFeatureExperiment),
[SpatialExperimentIO source](https://github.com/estellad/SpatialExperimentIO).

MultiAssaySpatialExperiment spatial layers/maps overlap representation needs,
not the refined sensitivity estimands. Inspected source 0.99.12 requires R≥4.6;
this host has R4.5.2, so only source review is claimed. Giotto/GiottoClass supplies
overlap/aggregation and multi-spatial-unit building blocks; the full R source
review did not identify the integrated A–H analysis.
[MultiAssaySpatialExperiment source](https://github.com/Genentech/MultiAssaySpatialExperiment),
[Giotto source](https://github.com/giotto-suite/Giotto).

These are source/API-fit judgments, not executable coercion certification. Actual
package adapters, CI, release compatibility and package check compliance belong
to later software development, not evidence blockers requiring API design to wait.
Extending a maintained package remains possible if maintainers prefer it, but
no such agreement exists and none is assumed.

## 15. Portability/adversarial tests

`144_validate_final_portability_cases.py` is a private independent scientific
oracle, not an implemented package API. Thirty-nine cases pass: numeric, very
large integer, character and UUID-like IDs; preservation of integer versus character
namespace distinctions; rejection of unsafe float/bool/missing IDs; empty censuses;
zero correspondence; exact identity/split/merge/complex; one/both-side dominance
ties; ambiguity; duplicate/unknown IDs; polygon overlap; holes; multipart/invalid/
empty geometry; boundary-inclusive touch; missing endpoint states; zero matched
mass/NA; feature identity/order mismatch; unit/origin/specimen mismatch; crop-edge
status; unchanged denominators under support loss. Production known-case checks
in both clean branches also pass.

The general feature contract fails on unexpected order or identity; the explicit
native-only fitted mask is not a license for silent query intersection. Coordinate
metadata disagreement fails rather than guessing a transform. Native IDs are
never spatially remapped. Invalid geometry remains explicitly ineligible.

Verdict **PASS WITH LIMITATIONS**: definitions/general edge cases are tested, but
this is not empirical validation on CosMx/MERSCOPE or executable Bioconductor S4
interoperability. Later implementation must carry these tests, not assume the
private oracle is a production library.

## 16. Final findings table

No confirmed CRITICAL defect in the qualified clean analyses. Findings below
include corrected implementation defects and enduring scientific qualifications.

| Severity | Exact source/logic | Problem and affected claim | Quantified effect | Fix / rerun / decision effect |
|---|---|---|---|---|
| MAJOR, resolved | `scripts/116_freeze_phase7_xenium_reference.R`, all-gene ≥3-cell assertion; strict relocated ROI 2 source | Fixed input panel mistaken for universally fit-able reference features; prevented replication | Slc6a3 zero in eligible native cells; strict fit aborts | User-approved native-only fixed-threshold mask in `151`; 247 fitted genes, all248 libraries/molecules unchanged. One amended reference/transfer fit, no segmentation rerun; no remaining API gate blocker |
| MODERATE, resolved additively | `scripts/118_compute_phase7_stability.py`, secondary_anchor_pair generic stability field | Pair agreement could be mistaken for incident-cell or mass-weighted agreement | ROI1 Baysor .914614 unweighted pair versus .944123 weighted pair-edge | `154` explicit pair table and specification separate all estimands; historical columns untouched; no numeric rerun needed |
| MODERATE, bounded | `116/117`, independently learned per-ROI states | 7 versus4 primary states can change numerical/rare effects; state IDs not common biology | ROI2 Baysor S/R .854467 local versus .881515 common reference; rare ProSeg merge signs change | Exact ROI1 frozen-reference transfer branch `148/150/152/153`, independently verified `155`; kept secondary. No universal topology penalty |
| MODERATE, bounded | `114` full components and `118` anchor-topology strata | Selected topology is not full-population census; weak edges join giant components | Rare ProSeg strata 1–7 cells; large selected gaps up to .247 | Primary all-edge/all-cell results; full component/undefined strata and selection fractions explicit. Reduced-stability claim only partially supported, not causal/global |
| MODERATE, bounded | `129/130` strict bbox-interior and local segmentation output | Geometric interior does not ensure independence from segmentation window | All1203 Baysor polygons inside ROI2, yet outside-molecule context unavailable | Specification distinguishes bbox status/window effects; no claim of full-tissue reconstruction. No prospective-outcome ROI changes |
| MODERATE, novelty qualification | `clustComp/R/flatVSflat.R`, intersection contingency graph | Shared-element weighted bipartite representation is not novel | Existing maintained generic functionality directly overlaps B | Restrict novelty to integrated mass/state/support/replicate/selection layer. NOVELTY CONDITIONAL; no priority claim |
| MINOR, corrected diagnostic | `147`, initial shortest-decimal float32 export to sf | Independent diagnostic rounded native coordinates when read as double | Initial ROI1 330/8112 Baysor hit counts changed; not a historical assignment defect | Preserve failed sample, exact float64/17-digit export; all4methods in bothROIs then agree exactly. Scientific pipeline unchanged |
| MINOR, corrected adapter metadata | `146`, sf/GDAL GeoJSON WGS84 default | Local µm coordinates wrongly tagged geographic by reader | Would invoke inappropriate geographic operations | Explicitly clear incorrect CRS without moving values; certified local-coordinate policy in specification |
| MINOR, resolved orchestration | `151/150/152`, relocated helper OUT and symlink producer lookup | Interrupted post-fit downstream execution/certification | Failed receipts preserved; no altered segmentation/model output | `153` resumes qualified paths and resolves producer hash by actual physical path. Completed/secondary PASS receipts authoritative; no resampling |
| NO ISSUE within tested scope | `141/142/146/155`, assignments/graphs/reference/stability/support | Tried to falsify mass, IDs, denominators, ties, projection and voting | Exact integer agreement; floating agreement within stated tolerance | Verified independently as described; does not certify absent provenance or general biology |

Initial failures are not counted as successful runs; completed qualification is
documented by `roi2_qualified/completed_execution_manifest.json` and the secondary
branch receipt. Historical outputs/scripts are unchanged. No silent cleanup or
replacement of adverse outcomes. Auxiliary Phase 7 calibration thread-environment
capture is less complete than final execution; this limits exact resource replay,
not scientific settings, and the new run records the complete environment/guard.

## 17. Final claim audit

| Claim | Adjudication | Bound |
|---|---|---|
| 1 Shared-molecule alternative-segmentation comparison | SUPPORTED | Certified same molecule universe, both clean ROIs |
| 2 Relationships lost by one-to-one geometric matching | SUPPORTED | Many-to-many mass information; not new graph mathematics or predictive superiority |
| 3 Directional denominators expose asymmetry | SUPPORTED | Independent assignment populations, not method accuracy |
| 4 Full versus selected-anchor stability can differ materially | SUPPORTED | Both signs possible; not universal optimism |
| 5 Structural disagreement does not uniformly imply state instability | SUPPORTED | Descriptive frozen expression states |
| 6 Some structural patterns associated with lower stability | PARTIALLY SUPPORTED | Baysor selected split/complex evidence; not global topology penalty |
| 7 Associations method/direction/ROI dependent | SUPPORTED | Rare/reference-dependent signs shown, one specimen |
| 8 Stochastic segmentation variability quantified explicitly | SUPPORTED | Three runs each ROI; conditional ARI with status accounting |
| 9 Useful sensitivity layer rather than truth score | SUPPORTED | Defined/descriptively demonstrated, no calibrated truth or clinical utility |
| 10 Defensible integrated Bioconductor software gap | PARTIALLY SUPPORTED | Bounded code review; established generic prior art; no complete A–H tool identified |

Machine-readable versions are `cross_roi/final_claim_audit.csv` and
`cross_roi/cross_roi_claims.csv`. Neither inconvenient cross-ROI results nor negative
anchor gaps were averaged away.

## 18. Remaining limitations

One specimen, two small ROIs, 248-gene targeted panel; no biological replication,
cross-platform empirical validation or general segmentation-quality ranking.
ROI boundaries truncate molecular context. ProSeg prior dependence and shared-
molecule dependence are not independent algorithmic evidence or causal prior tests.
Baysor threaded stochastic/executable behavior is not promised byte-identical on
rerun; observed outputs/provenance are frozen, while analytical/reference repeat
checks apply to the recorded environment. Alternative polygon PIP differs from
each algorithm's native internal assignment and native polygons are visualization
approximations. Ambiguity is measured, not resolved by assigning a nearest cell.
Rare topology effects and reference choices are unstable; no inferential p-values
or universal penalty. Winning-vote confidence is uncalibrated. Novelty review is
bounded, source assessed and non-exhaustive. The contract oracle is not a packaged
multi-platform library; current Bioconductor interoperability remains to implement.

These qualifications narrow future software claims but are not hidden open-ended
requests for another evidence phase. The failed mixed-generation breast geometry
does not re-enter validation. No historical result was repaired heuristically.

## 19. Final evidence gate

| Gate | Final classification |
|---|---|
| SCIENTIFIC REVALIDATION | PASS WITH LIMITATIONS |
| CROSS-ROI REPLICATION | PASS |
| CORE METHOD | VALIDATED |
| BIOLOGICAL SENSITIVITY | VALIDATED DESCRIPTIVELY |
| PORTABILITY TESTS | PASS WITH LIMITATIONS |
| NOVELTY | CONDITIONAL |
| BIOCONDUCTOR CONTRIBUTION TYPE | NEW PACKAGE |
| FINAL PACKAGE DECISION | GO TO API DESIGN |

Cross-ROI PASS means the central scoped sensitivity conclusions replicate, not
every candidate claim or numeric effect. The explicit partially replicated anchor
claim and rare-effect failures remain part of the evidence. Novelty CONDITIONAL
means software integration distinctness is defensible but mathematical/publication
priority is not proven; it is not a blocker for the narrowed design scope.

No unresolved prerequisite genuinely needs another execution/audit pass before
API design: the native-only technical feature-policy blocker was approved, tested
and documented in this pass. API design should start from the lossless molecule
universe, segmentation-specific assignment/status censuses, full weighted graph
and explicit mass/selection estimands, with optional frozen-reference propagation;
reuse existing experiment/geometry classes instead of inventing a spatial container.
Carry the scientific contract/adversarial cases into implementation. Normal later
software obligations—adapter/CI tests, package checks, licensing and release fit—
remain implementation work, not an additional pre-API scientific phase.

New scripts are **140–157**, listed with hashes in the final ledger: prospective
execution; independent graph/reference export and audit; maintained-source snapshots;
general portability oracle; sf checks; exact common-reference branch; adjacent
clustering source review; approved feature amendment and completion; cross-ROI
tables; additional numerical checks; novelty adjudication; final hash/gate closure.
All new scientific outputs are isolated under `results/final_pre_api/`: original
ROI2 generation/strict failure, qualified amended primary, secondary common-reference,
independent audits/geometry samples, code snapshots, novelty matrix, portability,
cross-ROI tables and closure manifests. The two root scientific documents are new.
No historical output or script 00–139 was overwritten and no package API work began.
