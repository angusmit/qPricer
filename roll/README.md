# roll/ — Roll discipline (Research OS R11)

## Purpose
Commodity backtests are most often wrong for one of two reasons: rolling on data the strategy could not have seen, or backtesting an abstract continuous series instead of the contracts you would actually trade. R11 makes roll discipline a **config + a single function**: which contract is **ACTIVE** as of a date, decided **DETERMINISTICALLY FROM AS-OF DATA ONLY** (through R9's door), with **roll events** recorded and the **continuous series demoted to an analytics-only derived view**.

**THE PRINCIPLE: trade the actual contracts the roll map names; never trade or signal off the continuous series.**

## Dependencies
A **LOW** layer: depends **ONLY on R9's door** (`.state.asof`) for as-of expiry / volume / prices; **independent of `curve/`** (placed after it for tidiness). Reads `state/` (and via it `data/`) downward only; opens nothing at import. **Additive — it does NOT rewire the existing readers.** Registered as an R2 `roll` capability (`rollEngine`, conforms `.contracts.verify`); carded via R5.

## Modules
- `roll.q` — `.roll.*`: the rule engine, roll events, and the analytics-only continuous view.

## Key API
- `.roll.active[asOf;comm]` → the active contract(s) as of `asOf`, from the door's as-of data ONLY: `active` (the primary `contractYM`), `weights` (`contractYM!weight` — one contract normally, two during a roll window), `reason`, plus `asOf`/`commodity`/`effectiveDate`/`ruleType`.
- `.roll.events[comm;fromDate;toDate]` → an `([] date; commodity; fromContract; toContract; reason)` table — one row each time the active contract changes over the window. Dates with no live contract (the data edge) are dropped.
- `.roll.continuous[comm]` → `([] date; contract; rawFront; adjusted; rolled)` — the **analytics-only** back-adjusted continuous series (see the warning below).
- `.roll.eventsPersist[hdbPath;events]` / `.roll.eventsOpen[hdbPath]` — persist / load the `rollEvents` splay (splayed, **unpartitioned**, `p#` on `commodity`, `.Q.en` + `set` — the same storage decision as the HDB; gitignored).

## Roll rules (`.cfg.rolls`, per-commodity with a `default` fallback)
- `days_before_expiry` **(default)** — hold the nearest-expiry contract whose `(expiry - asOf) > rollDays`; in the `W` days before that contract's own roll point, **blend** held→next so the roll is not a single-day jump.
- `volume_switch` — front vs next by **trailing as-of volume**; roll when the next contract's volume exceeds the front's.
- `fixed_calendar` — roll to next on/after the `rollDays`-th day of the front's expiry month.
- `oi_switch` — **DESIGNED but FLAGGED UNAVAILABLE.** The HDB carries **volume + expiry but NO open interest**, so requesting `oi_switch` errors rather than fabricating OI (the same honesty as R9's revision-awareness without fabricated revision columns).

## As-of-only by construction (the careful-reviewer point)
Every rule consumes ONLY `.state.asof[asOf;comm]`'s output (rows with `date<=asOf`), so a roll decision **cannot see future data**. The test plants a **future volume spike** (`date>asOf`) and asserts the past `volume_switch` decision is **unchanged** — point-in-time safety is inherited from the door, not re-implemented.

## The continuous series is ANALYTICS-ONLY (the third careful-reviewer point)
`.roll.continuous` back-adjusts at each roll using the **new contract's own move** (difference method) so the series has no gap — but back-adjustment **rewrites historical levels at every roll**, so the series is **NON-POINT-IN-TIME** and is for **charting / long-run vol ONLY**. Never trade or signal off it; the point-in-time-safe, tradable path is `.roll.active` (the actual contracts). R13 will flag any strategy that signals off this series. (The active-contract *choice* per date is still as-of-clean via `.roll.active`; only the price back-adjustment is retrospective.)

## Notes
- Config in `.cfg.rolls`. Carded (R5, `card_rollEngine`) with real failure modes: as-of-only decisions vs the temptation to peek; `oi_switch` unavailable on this data; the continuous series being mistaken for tradable.
- ADDITIVE / byte-identical: `roll/` is new and reads-via-the-door + computes; it edits no compute path. The existing suite stays byte-identical; +2 synthetic roll tests.
- Reserved-name discipline: parameter `asOf` (never `asof`), `comm` (never `commodity`); builtins (`next`/`prev`/`ratios`/`deltas`/`differ`/`first`/`last`/`bin`/`binr`/`get`/`key`/`value`) not shadowed.
- Demo: `apps/examples/roll_discipline.q` (real CRUDE — active contract, the roll-event cadence, continuous vs actual front, the as-of-only property).
