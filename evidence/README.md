# evidence/ — The Evidence Audit (Research OS R13)

## Purpose
The evidence-layer **gate** that closes the loop the door opened. R9's accessor **prevents** look-ahead by construction; R12's replay **records** the provenance (every step's as-of data access, the book, the fills, the rolls, the PnL components); R13 **verifies** it — a deterministic pass/fail audit that bites on a replay run record **before** the governance gates ever see its PnL. After R13, the only PnL the gates judge is PnL **proven** to correspond to something a desk could actually have traded.

## Dependencies
A layer **ABOVE `backtest/`** (it consumes R12's run record) and **BELOW `gov/`** (it reads the zone **config** `.cfg.gov.zones`, not gov's logic). Reads `state/` (the as-of universe), `backtest/` (the run-record shape), and `config/` (zones + tolerance) — all downward. Loads after `backtest/`, before `gov/`; opens nothing at import. **ADDITIVE: it changes no existing compute path, leaves `gov/` UNMODIFIED, and leaves `.workflow.run` (R7) UNCHANGED.** **INFRASTRUCTURE** (a gate, like the gov gates + the replay engine) — NOT registered as an R2 capability and NOT carded.

## Modules
- `evidence.q` — `.evidence.*`: the deterministic verifier + the fail-safe precondition wrapper.

## Key API
- `.evidence.audit[run]` → `` `pass`checks`reason `` — a pure, deterministic verifier of an R12 replay run record. `checks` is a dict checkName→boolean; **any** false ⇒ `pass=0b` with the `reason` naming the failed checks. The eight checks:
  - **lookAhead** — every step's `asOf` is its own date and no as-of slice saw data dated after `asOf` (`provDateTo<=asOf`); R9 prevents it, the audit verifies it.
  - **universe** — every fill's contract was live (expiry`>`asOf) as of its date.
  - **rollRespected** — the recorded `rollEvents` are exactly the transitions in the active-contract series (no fabricated/missing roll) and every roll is forward (no arbitrary backward roll).
  - **costsApplied** — a run with fills carries a positive transaction cost (unless `requireCosts` is off).
  - **fillsPresent** — no position change without a corresponding fill.
  - **bookTies** — the position book equals the cumulative fills.
  - **pnlTies** — `stepPnl ≈ positionPnl − totalCost − financingCost` per step, within `.cfg.evidence.reconcileTol` (the "is this PnL real" reconciliation).
  - **containment** — the run's dates ⊆ train+validate per `.cfg.gov.zones` (the holdout was not spanned — a light cross-check; the seal stays in gov, R3b).
- `.evidence.gatedRun[hypoId;run;runner;axis]` — the **fail-safe precondition wrapper**, mirroring `.cards.gatedRun`: audit the run; on **PASS** delegate to the gates (`.gov.runFull`); on **FAIL** return a reject verdict `evidenceFailed` (`tradeable=0b`) and **the gates NEVER run** (gov is not invoked). A run that cannot be **proven** clean is rejected, not judged.

## The audit BITES (the discipline)
An audit that always passes manufactures false confidence, so **every** check has a synthetic test that feeds a deliberately-contaminated run and confirms the audit catches it (`tests/evidence/test_audit_bites.q` — one contamination per check). Same discipline as R2 conformance / R5 `.cards.audit` / R6 stationarity.

## Notes
- Config in `.cfg.evidence` (`reconcileTol`, `requireCosts`). The containment boundary is computed from `.cfg.gov.zones`'s floor-cut formula + the HDB date axis (reading the **config**, not gov's logic; not a second holdout seal).
- Reserved-name discipline: `asOf` not `asof`, `comm` not `commodity`; builtins (`all`/`any`/`sum`/`sums`/`deltas`/`differ`/`min`/`max`/`first`/`last`/`value`/`within`) not shadowed. **q trap learned here:** qSQL `select`/`exec`/`update` on a **local** table inside a lambda throws `'assign` — the audit filters with plain vector ops (`steps where mask`) instead.
- Demo: `apps/examples/evidence_audit.q` (real CRUDE clean run PASS + synthetic contamination FAIL + `.evidence.gatedRun` rejecting the contaminated run while passing the clean one to the gates).
