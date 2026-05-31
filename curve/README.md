# curve/ — The curve engine (Research OS R10)

## Purpose
The dedicated curve capability: from the as-of data (obtained THROUGH R9's single door, so point-in-time is inherited), construct the clean curve + calendar spreads and compute the **rich derived features** — roll yield, slope, curvature, contango/backwardation classification — plus the **curve-shock operators** and **immutable daily snapshots**. It is the curve source the replay engine (R12), seasonality/carry (R14), PnL/risk attribution (R15) and the first real strategy (R16) all build on.

## Dependencies
A **LOW** layer: loads **after `state/`** (consumes `.state.asof` / `.state.build`) and **before** the consumers. Reads `state/` (and via it `data/`) downward only; opens no HDB/snapshots at import. **Additive — it does NOT rewire `regime/` or `.state.build`** (regime/ keeps its own slope/classification; rebasing it onto the curve engine is a later deferral). Registered as an R2 `curve` capability (`curveEngine`, conforms); carded via R5.

## Modules
- `curve.q` — `.curve.*`: the engine, the shock operators, and the immutable snapshot facility.

## Key API
- `.curve.build[asOf;comm]` → a curve object dict: `asOf`, `commodity`, `effectiveDate`, `curve` (clean, the `.data.hdb.curveAt` shape — `tenor`[years]/`price`/`contractYM`/`expiry`), `spreads` (adjacent calendar spreads), `features`:
  - `slope` — `(deferred-front)/front` (**matches regime/**), `classification` — `backwardation`/`contango`/`flat` by the **same flat threshold as regime/**, `curvature` — the butterfly second difference (`near - 2·mid + far`), `rollYield` — annualised `ln(Ffront/Fdeferred)` (positive in backwardation). Deterministic.
- `.curve.shock.parallel|slope|butterfly[curve;size]` (pure) — add a constant / a front-pivot tenor-proportional amount / a `4x(1-x)` belly bow. For the R15 scenario / bucketed-risk work.
- `.curve.snapshot[asOf;comm]` — compute via `.curve.build` and persist immutably to the in-memory `curveSnapshots`. `.curve.snapshotAt[date;comm]` / `.curve.snapshotHistory[comm]` read; `.curve.snapshotPersist`/`.snapshotOpen` set/get the splay.

## Agreement with regime/ (the precondition for a later rebase)
The slope and the backwardation/contango/flat classification use the **same convention + threshold** as regime/ — `.cfg.curve`'s `deferredIdx` and `flatThreshold` are **sourced from `.cfg.regime`** so they cannot drift — and `tests/curve/test_curve_engine.q` asserts `.curve.__features`'s classification equals `regime/`'s (`.regime.__labelPanel`) on the same curve. The two derivations agree; regime/ can be rebased onto this engine later with no surprise.

## Snapshots — splayed, UNPARTITIONED (matching the HDB decision)
`curveSnapshots` is stored the **same way as the HDB (§3)**: splayed, **unpartitioned**, `p#` on the `commodity` column, sym columns enumerated via `.Q.en`, written with `set` — **not date-partitioned** (~daily partition dirs would be slow + Windows-hostile). **Immutable:** a snapshot for a (commodity, date) is written once and never overwritten in place; recomputing must produce the identical curve (deterministic), so a re-snapshot is a no-op-if-matching and a **differing** re-write errors (catching non-determinism). The `curveSnapshots` dir is gitignored (like the HDB); tests use a synthetic table only.

## Notes
- Config in `.cfg.curve`. Carded (R5) with real failure modes: thin-liquidity tenors distorting slope/curvature; classification threshold sensitivity near the flat boundary; regime/ not yet rebased (two derivations coexist — mitigated by the agreement test); snapshot immutability assumes a deterministic curve.
- ADDITIVE / byte-identical: `curve/` is new and reads-via-the-door + computes; it edits no compute path; regime/ and `.state.build` are untouched. Demo: `apps/examples/curve_engine.q`.
- Reserved-name discipline: parameter `asOf` (never `asof`), `comm` (never `commodity`); builtins (`deltas`/`avg`/`dev`/`var`/`cor`/`cov`/`sum`/`first`/`last`/`differ`) not shadowed.
