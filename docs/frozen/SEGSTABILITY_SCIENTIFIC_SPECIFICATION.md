# SegStability consolidated scientific specification

Scope: descriptive sensitivity of alternative cell partitions of the same observed
molecule universe. This is a scientific contract, not a package API design. The
historical mixed-generation breast dataset is not validation evidence. No
segmentation, polygon, transferred state or native expression state is biological
truth. The final pre-API audit adjudicates empirical evidence separately.

## 1. Input and provenance requirements

One immutable evaluation universe U per specimen/ROI is required. Each retained
molecule has a unique, lossless original identity, feature identity, original
coordinates and filter metadata. Filtering is frozen before outcomes. A molecule
must occur at most once in U; duplicate identities fail. IDs are opaque, scoped by
specimen and segmentation; equal cell-ID spellings across segmentations do not
establish identity. Large integers must never pass through floating point. Explicit
lossless integer-to-decimal serialization is permitted with type and round-trip
evidence; heuristic relabeling or spatial ID repair is not.

Each segmentation has an explicit cell census, including zero-molecule and
zero-correspondence cells, and a molecule assignment status. Assigned means one
certified cell ID; unassigned and ambiguous/unresolved have no forced cell ID.
Ambiguous candidate sets/hit counts must be retained. Missing IDs, duplicate cell
IDs, assigned IDs outside the census and forced ambiguous assignments fail.
An empty segmentation is valid graph input, not a missing file or failed run.

Certified authoritative mask assignments may be used without polygon reassignment,
but their ID/coordinate linkage must be proven. If assignment uses geometry,
polygon/multipolygon rings, holes and object identity must be retained; invalid,
empty and non-positive-area objects remain in the census as explicitly ineligible.
No implicit repair, buffering, nearest-cell fallback, clipping or smoothing.
Boundary-inclusive point intersection is the clean validation policy: zero hits
unassigned, one hit assigned, multiple hits ambiguous. Interior hole points do not
intersect; hole-ring boundary points do. Invalid geometry is excluded from the
assignment index, never silently converted to valid geometry.

Common specimen, units, origin, axis direction and every transformation must be
declared and verified before spatial operations. Different CRS/units or unverified
registration fail; coincident numerical ranges do not certify alignment. GeoJSON
readers can automatically attach WGS84 to local micrometre coordinates. A certified
local-coordinate adapter must explicitly remove that incorrect metadata without
transforming coordinate values. It must not treat micrometres as longitude/latitude.

Source URL/version/retrieval date, input/software/source hashes, command, model,
parameters, thread/runtime environment, observation window and output hashes are
required. Unknown provenance must remain unknown; missing required evidence closes
the corresponding gate. Resource limits may reduce a prospective spatial window,
not retrospectively tune scientific settings for favorable results.

## 2. Weighted bipartite correspondence graph

For two segmentations X and A, nodes are their full explicit censuses, separately
namespaced. For i in X and j in A, edge weight w(i,j) is the integer number of
molecules in U uniquely assigned to BOTH i and j. Only positive edges are stored.
Ambiguous/unassigned molecules are not forced onto edges. Geometry is not needed
when assignment identities are already certified.

Independent denominators N_X(i) and N_A(j) count ALL uniquely assigned molecules
on the respective side, including molecules unmatched or unresolved on the other
side. They are not row/column sums of the shared-edge contingency. Directional
fractions are w/N_X and w/N_A. Each fraction lies in [0,1]. Cell matched mass M_i
is incident-edge sum; unmatched mass is N_i-M_i. Zero-edge cells remain explicit.
Zero own denominator is non-estimable, not zero or perfect agreement. Global
assigned+unassigned+ambiguous counts equal |U| separately per segmentation. Global
ambiguity is not conflated with unmatched assigned mass or missing endpoint labels.
Cross-assignment status tables are needed to explain these overlapping conditions.

## 3. Primary structural summaries

Partner count is positive-edge degree. Dominance maximizes integer w on the
COMPLETE observed graph; for a given cell, its independent denominator is constant,
so maximizing its directional fraction is equivalent. All exact maximum-weight
ties are retained. No edges means no dominant partner and zero tie count; one
maximum means unique dominance; more than one means unresolved dominance.
Reciprocal dominance requires a unique maximum on both endpoints of the SAME
edge. A counterpart tie is not reciprocal. No arbitrary ordering tie-break.

Connected components include isolated census nodes. Component cardinalities are
reported as (n_X,n_A): (1,1) one-to-one; (1,>1) split; (>1,1) merge; (>1,>1)
many-to-many; a one-sided isolated node unmatched. These are component summaries,
not the secondary reciprocal-anchor classes. A single weak edge can connect a
giant component; component labels are not a census of biological split/merge
events. All-cell directional matched fractions, partner distributions, ties,
unmatched and ambiguous mass are primary. Complete graph summaries precede selection.

## 4. Secondary reciprocal-anchor topology

The clean analysis eligibility policy is valid positive-area geometry, strict
bbox interior to the frozen inclusive rectangular ROI on BOTH endpoints, and
unique reciprocal dominance computed on the COMPLETE graph. Full source polygons
are retained, not clipped. Touching the ROI boundary is crop-edge, not strict
interior. Missing geometry yields unknown status, never inferred from molecule
locations. Locally generated polygons contained in the ROI do not certify that
segmentation would be unchanged if outside-ROI molecules were available; window
effects remain a limitation distinct from bbox status.

For an eligible reciprocal edge (i,j), fan-out(i) is the number of eligible
interior A cells whose unique complete-graph dominant X partner is i; fan-in(j)
is the number of eligible interior X cells whose unique complete-graph dominant
A partner is j. Both endpoints of each contributing dominance relation must be
interior. Classify the anchor edge as:

| fan-out | fan-in | secondary class |
|---|---|---|
| 1 | 1 | one_to_one_like |
| >1 | 1 | split_like |
| 1 | >1 | merge_like |
| >1 | >1 | complex |

These classes describe dominance relationships, not verified error or causal
biology. No reciprocal eligible anchor means topology undefined. Non-reciprocal
matched cells, tied dominance, reciprocal non-interior cells and zero-correspondence
cells retain explicit strata; never force them into these four classes. Support
thresholds do not silently redefine primary dominance or anchor eligibility.

## 5. Frozen expression-state endpoint

Fit exactly ONE reference using native Xenium expression only for each ROI in the
clean replication. Alternative expression must not affect eligibility, normalization,
scaling, PCA, graph, states or parameter selection. Counts use ONLY retained ROI
molecules, not whole-cell counts including molecules outside the window.

Clean empirical recipe: retain all 248 Gene Expression input features in accepted
matrix order and all evaluation molecules. Native cells require at least five
detected genes. Before reference fitting, apply the fixed >=3-eligible-native-cell
feature rule and freeze the resulting fitted gene order for every alternative in
that ROI. ROI 1 fits all 248 genes; ROI 2 fits 247, excluding absent Slc6a3. This
explicit native-only feature-policy amendment was approved in the final pass;
the original all-features assertion failed and is preserved, not hidden.
Normalize fitted genes as log2(1+5000*count/library), where library includes ALL
248 input features, including any not fitted. Use native training
feature means/sample SD, no cell z-scoring; full eligible-native PCA with 30 saved
dimensions, PCs 1–15 used. Record full input order, fitted order and every feature
eligibility count. Non-finite/zero SD among fitted features or insufficient eligible
cells/features fails; there is no silent or alternative-informed feature selection.
For other datasets, the input feature population and native-only technical policy
must be prospectively specified, not selected using stability outcomes.

Reference graph: exact k=30 non-self nearest neighbors by identity, candidate
union of directed kNN edges, weight=shared-neighbor count/30; no diagonal edges.
Weighted undirected modularity Leiden, resolution 0.15, 100 iterations, seed 2026
and frozen package/RNG versions. Sensitivity resolutions 0.10/0.30 use the same
PCA/graph; repeated-seed memberships are checked. Freeze/hashes record genes,
all/eligible cells, transforms, means/SDs, loadings/scores/PCs, graph/kNN, states,
resolution/seed, transfer parameters, versions, input hashes and dependencies.

Query alternative counts must have identical full input feature identity/order;
the explicitly frozen fitted-feature mask must produce identical fitted order.
At least five detected FITTED features are required. Full-panel library sizes
remain unchanged. Normalize and project using only frozen
reference statistics/loadings. Exact reference-neighbor votes, k=30 primary,
k=15/50 sensitivity. A tied winning vote OR distance tie between kth and (k+1)th
neighbor within 1e-12*max(1,d_k) stays unresolved. No forced prediction, post hoc
confidence threshold or alternate-reference fitting. Record top/second vote,
margin, nearest-reference distance, status and unresolved cell/molecule counts.
Winning-vote fraction is a diagnostic, not a calibrated probability. Reference
states are Xenium-reference expression states, never true cell types. Independently
learned ROI reference-state integers are not aligned biological identities across
ROIs; cross-ROI comparisons assess estimand/qualitative robustness, not label identity.

## 6. Primary biological estimands and missing-state accounting

For each edge, labels are resolved only if BOTH endpoint states are available.
Let M=sum(w), R=sum(w for resolved edges), S=sum(w for resolved same-state edges).
Report S/R (conditional same-state agreement), S/M (same-state fraction of all
matched mass), S/sum(N_X), S/sum(N_A), M-R unresolved-label mass, each side's
sum(N)-M unmatched assigned mass, and global ambiguous assignment mass. Missing
labels are neither disagreements nor matches. All denominators are explicit;
zero denominators return non-estimable/NA.

For each cell, M_i/R_i/S_i sum ALL incident edges; primary directional total-mass
stability is S_i/N_i and conditional stability is S_i/R_i where estimable. Report
mass-weighted totals AND equal-cell means over their explicitly estimable
populations. Zero-molecule, zero-correspondence, endpoint-unavailable and
non-estimable cells are counted, not dropped. Strata use the same definitions
and original denominators. Neither matched mass nor resolved matched mass is a
substitute for the full eligible census.

## 7. Secondary anchor estimands and selection accounting

Keep three distinct quantities: (1) unweighted same-state agreement over
resolved reciprocal-anchor PAIRS; (2) Xenium-anchor CELL-population stability
using all incident edges of selected native cells; (3) alternative-anchor
CELL-population stability using all incident edges of selected alternative cells.
Additionally record selected-anchor EDGE mass and weighted pair-edge S/R, explicitly
different from (1). Unweighted pair agreement must not be named molecule-weighted
conditional stability.

Report anchor share of all census cells, incident matched mass and selected
pair-edge mass separately, anchor/non-anchor matched stability and signed gaps.
An optimistic anchor result is possible, not guaranteed. Population-flow conditions
must identify whether they are nested filters or non-nested census predicates;
retention percentages cannot imply a sequential exclusion mechanism if predicates
are non-nested. Endpoint missingness, ties, crop edges and zero matches remain visible.

## 8. Support sensitivity

Prespecified directional thresholds are 0%,1%,2.5%,5%,10%,20%. For side X retain
edges with w/N_X>=threshold; for side A use w/N_A. These are separate estimands,
not an OR/AND union unless a different policy is explicitly prospectively named.
Assess full observed support, both-interior support and selected anchor-pair edges.
Preserve complete-graph own denominators, dominance and node census. Report retained
edges, positive-degree retained cells, retained shared mass/fraction, excluded
matched mass, unresolved label mass, S/R and S/N together. A higher conditional
agreement bought by mass loss is not evidence of improved whole-population stability.

## 9. Dependence, stochastic runs and geometric comparator

All comparisons are paired transformations of the same molecules/tissue.
ProSeg native-prior/nucleus initialization is dependence, not independent replication.
Baysor no-native-prior does not eliminate shared-data dependence. Three stochastic
ProSeg runs are separate units of computational repeatability; cell IDs must never
be matched literally across runs. Match transcript identities, report status
contingencies, jointly uniquely assigned population and label-invariant partition
agreement, graph/state/support/rare-topology variability separately. Do not average
away inconsistent signs. Differences between prior-dependent and no-prior methods
do not isolate the causal effect of using a prior.

Use certified local coordinates for IoU/Jaccard (same geometric quantity), Dice,
directional area overlap and maximum-total-IoU one-to-one matching. Report all
positive geometry pairs as well as molecule edges, geometrically matched pairs
with no shared molecules and molecule edges omitted by one-to-one matching.
Native visualization polygons are not exact masks. Distinct retained information
does not prove predictive superiority or a new segmentation-quality truth score.

## 10. Inference boundary

Results are descriptive segmentation/reference-state sensitivity in observed ROIs.
Cells/molecules are dependent observational units, not biological replicates.
The biological inference unit is specimen/sample; disjoint ROIs of one specimen
test spatial robustness, not independent biological generalization. No cell-level
population p-values, universal topology penalty, causal topology effect, true cell
type assertion, segmentation-accuracy ranking or universal novelty claim is justified.
Portability contract tests certify definitions/edge cases, not empirical success
on every platform. A later software implementation must preserve this contract
while reusing existing Bioconductor infrastructure; representation alone is not
the claimed analytical contribution.
