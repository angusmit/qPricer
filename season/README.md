# season/ — Curve/spread seasonality (Research OS R14)

## Purpose
An evidence-layer **feature** capability: the curve/spread **seasonality** the first real strategy (R16) is built on. The seasonally-adjusted **same-calendar-month z-score** of the front-deferred spread IS the signal for R16's calendar-spread mean-reversion. Computed **causally** through R9's door, so the seasonal statistics are point-in-time by construction.

## Dependencies
A **LOW** evidence-tier layer: loads **after `curve/`** (reuses R10's front/deferred/slope **convention** + `deferredIdx`) and reads `state/` (R9) downward; opens nothing at import. **ADDITIVE** — reads + computes, edits nothing; does **not** rewire `.state.build`, and is **distinct from** the existing `signals/` seasonality (`.commodity.seasonality.*`, the alpha-library higher in the stack — no collision). Registered as an R2 `season` capability (`curveSeasonality`, conforms); carded (R5).

## Modules
- `season.q` — `.season.*`: the causal seasonal feature engine.

## Key API
- `.season.features[asOf;comm]` → a feature dict: `sameMonthZ` (the z of the current front-deferred spread vs its **same-calendar-month** history up to `asOf` — the R16 signal), `sameContractMonthZ` (vs the same front delivery-month historically), `seasonalFactor` (the calendar-month mean front level, causal), `deseasonalised` (front level ÷ seasonal factor), `seasonalSlope` (the slope minus its same-month mean), plus `spread`/`frontLevel`/`slope`/`nSameMonth`/`lowConfidence`.
- `.season.__series[asOf;comm;deferredIdx]` — the causal per-date front-deferred spread/slope series (one pass over the as-of slice).

## The point-in-time discipline (the key reviewer point)
The same-month z on date D uses **ONLY same-calendar-month observations up to D** (trailing/causal), obtained through R9's door (`.state.asof` returns `date<=asOf`). A **full-sample** seasonal statistic is look-ahead (it standardises a past observation by future same-months) — the silent killer of seasonal backtests. The test plants a **future** same-month observation and asserts a past z is **unchanged**. Below a minimum same-month N (`.cfg.season.minObs`) the z is `null` + `lowConfidence` rather than a spurious value (thin same-month history early in the data).

## Notes
- Config in `.cfg.season` (`deferredIdx` = the seasonal spread's deferred leg, default 2 = M1-M3; `minObs`). Carded with real failure modes (thin same-month history early → noisy z, mitigated by min-N; seasonal patterns break in structural shifts; causal trailing stats use less data than full-sample).
- Reserved-name discipline: `asOf` not `asof`, `comm` not `commodity`; calendar month via `(\`month$date) mod 12` (the `.mm` accessor is **atom-only** — it errors on a date vector); qSQL `select`/`exec` on a local table inside a lambda throws `'assign` (skill #24) — the series is built with plain mask+`where`+`group`.
- Demo: `apps/examples/season_carry.q` (real CRUDE: the causal same-month z vs the full-sample look-ahead contrast).
