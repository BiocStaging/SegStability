# SegStability Codex Instructions

## Project mission

SegStability is a Bioconductor package for **descriptive segmentation-sensitivity analysis** in molecule-resolved spatial transcriptomics.

It compares alternative cell segmentations of the **same immutable molecule universe**, quantifies their structural correspondence, and measures how segmentation differences propagate to downstream categorical cell-state conclusions.

It is **not**:

- a segmentation algorithm
- a segmentation-quality truth score
- a cell-type caller
- a method for deciding which segmentation is best
- a biological-ground-truth system
- a Baysor/ProSeg wrapper
- a replacement spatial-data container

Scientific correctness, lossless identity handling, reproducibility, and auditability take priority over convenience.

---

## Authority order

When instructions conflict, follow this order:

1. `docs/frozen/SEGSTABILITY_SCIENTIFIC_SPECIFICATION.md`
2. `docs/frozen/SEGSTABILITY_API_SPECIFICATION.md`
3. `docs/frozen/SEGSTABILITY_IMPLEMENTATION_PLAN.md`
4. `docs/frozen/SEGSTABILITY_FINAL_PRE_API_AUDIT.md`
5. this `AGENTS.md`
6. applicable installed skills
7. general coding conventions

The API and implementation plan are **FROZEN FOR v0.99 IMPLEMENTATION**.

Do not redesign the architecture or reopen scientific validation unless implementation exposes a concrete contradiction. Document that contradiction before proposing a change.

---

## Working method

Before substantive work:

1. read the relevant frozen specification sections;
2. identify the current implementation phase and prerequisite gate;
3. inspect existing code/tests before modifying anything;
4. make a short implementation plan;
5. implement the smallest scientifically correct increment;
6. run the relevant tests/checks;
7. verify the result before claiming completion.

Do not jump ahead to later modules.

Follow the implementation sequence in `SEGSTABILITY_IMPLEMENTATION_PLAN.md`.

---

## Skills and tools

Use installed skills/plugins when they materially improve correctness or efficiency.

Especially prefer:

- `bioc-pkg-dev`
- `bioc-pkg-review`
- `analyze-r-package`
- `developing-r-packages`
- `writing-r-code`
- `test-driven-development`
- `systematic-debugging`
- `verification-before-completion`
- `scientific-critical-thinking`
- `bioc-howto`

Use scientific/bioinformatics skills when relevant.

Use `caveman`, `ponytail`, or similar context/token-management skills when they help reduce unnecessary context or repeated work.

Do **not** invoke skills mechanically. Use only what helps the current task.

Do not create subagents unless:
- the task clearly benefits from independent parallel work, and
- the expected improvement justifies the extra context/token cost.

Prefer one focused agent working carefully over unnecessary orchestration.

---

## Token and context discipline

Usage efficiency matters.

- Read only the files/sections needed for the current task.
- Avoid repeatedly loading entire long specifications when targeted sections suffice.
- Do not regenerate unchanged documents or code.
- Do not produce long explanations unless needed for a decision or audit.
- Keep implementation notes concise.
- Prefer targeted searches over broad repository scans.
- Avoid duplicate research already settled in frozen documents.
- Do not spawn subagents for routine coding, testing, or documentation.
- Reuse established project decisions instead of reconsidering them.

Spend tokens where they improve scientific or implementation correctness.

---

## Scientific rules that must never be violated

- One immutable molecule universe per analysis context.
- Molecule IDs are opaque and must remain lossless.
- Never pass scientific identifiers through ordinary floating-point representation.
- Never infer molecule identity from row order.
- Never heuristically repair IDs.
- Assigned, unassigned, and ambiguous are distinct states.
- Never force ambiguous molecules to the nearest cell.
- Explicit cell censuses include zero-molecule and zero-correspondence cells.
- Independent denominators come from each segmentation's complete uniquely assigned population.
- Never derive those denominators from correspondence-edge row/column sums.
- Dominance is calculated on the complete graph using integer shared-molecule edge mass.
- Preserve all exact dominance ties.
- Connected-component classes and reciprocal-anchor topology are different concepts.
- Support filtering must never redefine the primary graph, original denominators, dominance, components, or census.
- Missing biological states are unresolved, not disagreements.
- Stochastic segmentation runs are aligned through molecule identity, never same-valued run-local cell IDs.
- Geometry must never be silently repaired, buffered, clipped, transformed, or otherwise altered.
- Unknown provenance remains unknown.
- Results are descriptive sensitivity analyses, not biological truth or segmentation-accuracy rankings.

If convenience conflicts with these rules, correctness wins.

---

## Resource constraints

Primary development machine:

- Lenovo Yoga C940
- 16 GB nominal RAM
- historically about 12.7 GiB Linux/WSL physical RAM available to the environment
- NVIDIA GTX 1650, 4 GB VRAM
- limited-memory laptop development environment

Be conservative with compute.

Before large jobs:

1. inspect available **physical RAM**;
2. do not count swap as usable working memory;
3. estimate expected memory use;
4. set a conservative ceiling;
5. monitor the process;
6. stop or skip if unsafe.

Do not run expensive jobs simply because they are possible.

Avoid:

- dense cell × cell matrices
- unnecessary full-data copies
- unnecessary geometry loading
- million-scale tests when tiny deterministic fixtures answer the question
- GPU work unless scientifically or computationally justified

The 5M-molecule profiling target in the implementation plan is optional on constrained hardware.

---

## Testing discipline

Use test-driven development for scientific invariants.

Port the historical adversarial oracle incrementally alongside implementation.

For each phase:

- write/port relevant tests first or concurrently;
- implement the minimum required behavior;
- run focused tests;
- then run broader package checks at the defined milestone.

Never weaken a scientific test just to make implementation pass.

Expected major checks include:

```text
R CMD build
R CMD check
R CMD check --as-cran
BiocCheck
```

Interpret `--as-cran` notes intelligently; package correctness is ultimately judged in the current Bioconductor devel environment with `R CMD check` and `BiocCheck`.

Do not claim completion without verification.

---

## Dependency discipline

Keep the core lightweight.

Do not add a dependency because it is convenient or commonly used.

Every dependency must have a concrete implementation need.

Core assignment-based analysis must work without `sf`.

Do not introduce platform-specific readers that duplicate existing Bioconductor import infrastructure.

Reuse established Bioconductor classes and generics whenever appropriate.

---

## Coding discipline

- Follow Bioconductor/R conventions.
- Prefer small internal helpers with clear invariants.
- Keep the public API exactly within the frozen specification.
- Do not export arithmetic helpers unnecessarily.
- Do not expose internal slots as user API.
- Avoid premature optimization.
- Optimize only after correctness is established and profiling identifies a real bottleneck.
- Preserve deterministic ordering where specified.
- Fail loudly on scientifically meaningful invalid input.

Comments should explain scientific or architectural reasons, not obvious syntax.

---

## Authorship and attribution

Repository and package authorship must reflect the actual human contributors.

Use:

- Rama Krishna Pudota as the author/maintainer where appropriate

Do not add:

- Codex attribution
- ChatGPT attribution
- AI-generated-by comments
- AI co-author trailers
- assistant names in AUTHORS, DESCRIPTION, commits, documentation, source headers, or package metadata
- signatures such as "Generated by Codex" or "Created with ChatGPT"

AI tools are implementation aids, not package authors or maintainers.

---

## Local implementation records

For every substantive implementation work package, maintain a concise local audit
record under:

`docs/implementation/`

These records should capture actual changes, tests, checks, decisions, deviations,
gate outcomes, and the next allowed implementation boundary.

`docs/implementation/` is intentionally local-only and excluded through
`.git/info/exclude`.

Do not commit or push implementation audit logs.

Keep `CURRENT_STATUS.md` concise and local. Permanent repository documentation
should contain only information that belongs in the package/project record.

Do not rely on chat/session history as the only record of implementation work.

---

## Scope control

Do only the current requested phase.

Do not opportunistically implement:

- later API modules
- new plots
- extra adapters
- continuous biological endpoints
- alternate geometry models
- disk backends
- Arrow/HDF5 support
- new S4 classes
- new scientific metrics

unless the implementation plan has reached that phase or a concrete blocker requires them.

If a useful idea is out of scope, note it briefly rather than implementing it.

---

## Completion standard

A task is complete only when:

- it matches the frozen scientific/API contract;
- relevant tests pass;
- no known invariant is violated;
- no unnecessary dependency or computation was introduced;
- results have been verified;
- the next implementation gate remains clear.

When uncertain, prefer a small correct step over a large speculative one.
