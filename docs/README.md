# SegStability Documentation Map

This directory contains the authoritative project documents used during SegStability v0.99 implementation.

## `frozen/`

### `SEGSTABILITY_SCIENTIFIC_SPECIFICATION.md`
Defines the scientific contract.

Use it for:
- molecule-universe rules
- assignment semantics
- correspondence definitions
- denominators
- dominance/ties/components
- biological-state estimands
- support sensitivity
- stochastic-run interpretation
- geometry rules
- inference boundaries

This document defines **what must remain scientifically true**.

---

### `SEGSTABILITY_API_SPECIFICATION.md`
Defines the frozen software architecture and public API for v0.99.

Use it for:
- S4 class responsibilities
- slot schemas
- public functions
- accessors
- identifier representation
- fingerprint model
- geometry ownership
- dependencies
- interoperability
- memory/scaling strategy

This document defines **what the package must look like**.

Status: **FROZEN FOR v0.99 IMPLEMENTATION**

---

### `SEGSTABILITY_IMPLEMENTATION_PLAN.md`
Defines the implementation order and phase gates.

Use it for:
- current implementation phase
- prerequisite gates
- test-porting sequence
- package-check milestones
- resource constraints
- when later modules may begin

This document defines **how the package should be built**.

Status: **FROZEN FOR v0.99 IMPLEMENTATION**

---

### `SEGSTABILITY_FINAL_PRE_API_AUDIT.md`
Records the final scientific validation and evidence review that justified moving to implementation.

Use it for:
- empirical-validation context
- known limitations
- novelty positioning
- prior-art context
- portability/adversarial evidence
- interpretation boundaries

This is primarily an **evidence/audit reference**, not an implementation specification.

---

## Authority order

When documents appear to conflict, use:

1. `SEGSTABILITY_SCIENTIFIC_SPECIFICATION.md`
2. `SEGSTABILITY_API_SPECIFICATION.md`
3. `SEGSTABILITY_IMPLEMENTATION_PLAN.md`
4. `SEGSTABILITY_FINAL_PRE_API_AUDIT.md`

Project-level instructions in `../AGENTS.md` govern how Codex should work with these documents.

Do not modify files under `docs/frozen/` during normal implementation.

If implementation exposes a genuine contradiction, document the exact issue first rather than silently changing a frozen document.
