# SegStability — Codex Start Here

This repository is implementing the frozen SegStability v0.99 architecture.

Before making changes:

1. Read `AGENTS.md`.
2. Read only the relevant sections of:
   - `docs/frozen/SEGSTABILITY_API_SPECIFICATION.md`
   - `docs/frozen/SEGSTABILITY_IMPLEMENTATION_PLAN.md`
3. Consult:
   - `docs/frozen/SEGSTABILITY_SCIENTIFIC_SPECIFICATION.md`
   - `docs/frozen/SEGSTABILITY_FINAL_PRE_API_AUDIT.md`
   only when the current task needs scientific-definition or evidence context.

Do not redesign the API or reopen scientific validation unless implementation exposes a concrete contradiction.

## Current implementation boundary

Implementation must proceed in this order:

1. development-environment gate
2. authoritative package-name collision check
3. repository/package skeleton
4. lossless-ID exactness laboratory
5. universe/fingerprint utilities
6. `SegmentationPartition`
7. `SegmentationSet`
8. `SegmentationComparison`
9. `compareSegmentations()`
10. remaining modules according to the frozen implementation plan

Do not jump ahead.

## Immediate task

At the beginning of a fresh Codex session:

1. inspect repository state;
2. identify the current unfinished gate;
3. summarize what that gate requires;
4. make a short plan;
5. execute only that gate or requested subtask;
6. run the relevant verification before claiming completion.

## Important constraints

- preserve molecule identity exactly;
- never pass scientific IDs through ordinary double;
- never infer identity from row order;
- do not weaken tests to make code pass;
- do not create dense cell-by-cell correspondence matrices;
- do not make `sf` mandatory for core workflows;
- do not run unnecessary large jobs on the development laptop;
- do not introduce new public API beyond the frozen specification;
- do not spawn subagents unless clearly justified.

Use installed Bioconductor, scientific, testing, debugging, and token-management skills when they materially improve the current task.

Keep work focused, scientifically defensible, and resource-aware.

## Source-of-truth order

If anything conflicts:

1. scientific specification
2. API specification
3. implementation plan
4. final pre-API audit
5. `AGENTS.md`
6. installed skills
7. general coding preferences

When uncertain, stop at the smallest concrete contradiction rather than improvising around the frozen contract.
