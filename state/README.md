# state/ — Point-in-time market state + as-of discipline (Research OS R9)

## Purpose
The **first milestone of the evidence layer** and the chokepoint the rest of it stands on: **one as-of accessor that is the only door to history**, plus a first-class **Market State** object every later foundation milestone consumes. R1–R8 built a rigorous *judge*; the evidence layer is the realistic evidence the judge is handed — and *a rigorous judge is only as honest as the evidence*.

## Why one door
Point-in-time discipline fails when it is a rule every future module must remember. Make it **one function to audit** instead. The accessor **prevents** look-ahead by construction (it returns only rows knowable as of `asOf`) and logs **provenance**; R13's evidence audit later **verifies** by reading that provenance. Get the door right and the evidence layer inherits point-in-time safety for free.

## Dependencies
A **LOW** layer: loads **after `data/`** (it WRAPS the unchanged `.data.hdb.*` query layer — `curveAt`/`curveHistory`/`dates`) and **before** the consumers (`regime/` etc.). Reads `data/` downward only; opens no HDB at import (the accessor requires the HDB already opened, like the query layer). **Additive — it does NOT rewire the existing direct readers** (`regime/`, `backtest/` keep their current HDB access; rebasing them onto the accessor is a later deferral). Registered as an R2 `state` capability (`marketState`, conforms); carded via R5.

## Modules
- `state.q` — `.state.*`: the as-of accessor, the Market State builder, the provenance log, and the no-look-ahead invariants.

## Key API
- `.state.asof[asOf;comm]` → `` `data`provenance `` — the **only door to history**. Returns ONLY rows knowable as of `asOf` (`date<=asOf`) via the existing query layer, and appends a provenance stamp to the append-only `.state.provLog`. **Revision-aware by design**: if the source carries a `validDate` column it takes the latest revision with `validDate<=asOf` per observation date; futures carry no revision columns, so they take the simple `date<=asOf` path (we do **not** fabricate revisable data we lack). Parameter is **`asOf`, never `asof`** (the q as-of-join keyword).
- `.state.build[asOf;comm]` → the **Market State** dict: `asOf`, `commodity`, `effectiveDate` (latest trading date ≤ asOf), `curve` (the as-of curve, the same shape `.data.hdb.curveAt` returns — consistent with what `regime/` sees), `spreads` (adjacent calendar spreads), `features` (an empty placeholder the curve engine R10 + season/carry R14 fill — declared now so the contract is stable), `universe` (contracts NOT expired as of asOf), `refs` (cost-model + roll-rule handles; the real roll map is R11), `provenance`. Assembled **lazily** per (asOf, commodity).
- `.state.invariant.*` — the no-look-ahead invariants the tests assert and R13 will reuse: `asofRespected` (no row `date>asOf`), `universeLive` (the as-of curve carries no expired contract), `provenanceConsistent` (the stamp matches what was returned).

## Notes
- Config in `.cfg.state` (universe/liquidity filters: `expiryStrictlyAfter`, `minVolume`). Carded (R5) with real failure modes: revisable data not yet handled (only futures' simple path is exercised); the existing direct readers not yet rebased; the provenance log is advisory until R13 consumes it; the universe filter assumes the HDB expiry field is correct.
- ADDITIVE / byte-identical: `state/` is new and wraps/reads; it edits no compute path. Demo: `apps/examples/market_state_asof.q`.
- **Unblocks R10–R16**: the curve engine builds the state, roll trades through it, the replay loop folds over it, seasonality/carry fill its features, the first real strategy is judged on it; R13 verifies the provenance the accessor logs — the door prevents look-ahead, the audit proves it.
