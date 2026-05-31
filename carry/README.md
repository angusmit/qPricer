# carry/ — Carry / storage economics (Research OS R14)

## Purpose
An evidence-layer **feature** capability: the **carry / storage economics** the carry strategies are built on — implied carry, convenience yield, cash-and-carry fair value, the carry signal, and an inventory-tightness proxy. Derived from R10's curve (through R9's door), so it is point-in-time by construction.

## Dependencies
A **LOW** evidence-tier layer: loads **after `curve/`** (reuses `.curve.build`'s curve + annualised `rollYield`) and reads `state/` via it, downward; opens nothing at import. **ADDITIVE** — reads + computes, edits nothing. Registered as an R2 `carry` capability (`carryEconomics`, conforms); carded (R5).

## Modules
- `carry.q` — `.carry.*`: the carry/storage feature engine.

## Key API
- `.carry.features[asOf;comm]` → a feature dict:
  - `impliedCarry` — the curve-implied carry/roll, annualised `ln(Ffront/Fdeferred)/T` (**positive in backwardation**, negative in contango; reuses R10's `rollYield`).
  - `convenienceYield` — `(r+storage) − impliedCarry` (**assumption-dependent**).
  - `fairValueCarry` — the cash-and-carry no-arb deferred price `Ffront·exp((r+storage)·T)` (compare to `actualDeferred`; **assumption-dependent**).
  - `carrySignal` — `impliedCarry − (r+storage)`: the roll yield net of the assumed cost of carry (the edge signal; flips sign with the curve).
  - `inventoryTightness` — a **proxy** from the degree of backwardation (no real inventory data).
  - plus `front`/`deferred`/`tenorGap`/`classification`/`costOfCarry`/`assumptionDependent`.

## Honest data scope (the reviewer point)
The implied carry comes from the **curve** (we have it). The convenience yield + cash-and-carry fair value need a rate `r` + a storage cost, which we **supply from `.cfg.carry`** (config, **not** market-observed) — so those outputs are **assumption-dependent** and listed in `assumptionDependent`. The inventory-tightness is a **proxy** from the degree of backwardation; we have **no** real inventory data, so it is labelled a proxy and inventory is never fabricated (the same discipline as R11's no-open-interest).

## Notes
- Config in `.cfg.carry` (`riskFreeRate`, `storageCost`); carry reuses `.cfg.curve.deferredIdx` (the CL1-CLk leg) so `impliedCarry == rollYield`. Carded with real failure modes (convenience yield / fair value depend on the assumed r+storage; inventory-tightness a proxy; cash-and-carry assumes storability).
- Reserved-name discipline: `asOf` not `asof`, `comm` not `commodity`; builtins (`exp`/`log`/`sum`/`avg`) not shadowed; mind right-associativity (`front*exp costOfCarry*dt` = `front*(exp(costOfCarry*dt))`).
- Demo: `apps/examples/season_carry.q` (real CRUDE: the carry signal flips sign — backwardation vs contango).
