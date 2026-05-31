# regime/ — Market State Engine (ARCHITECTURE.md Part II §11.1-11.2, R1)

## Purpose
Deterministic **measurement** of a per-(commodity,date) regime fingerprint from the existing HDB futures data — curve state, realized-vol state, liquidity state, roll proximity, seasonal phase — so any backtest result can be broken down **by regime** (the difference between "0.8 Sharpe" and "+1.4 in backwardation, negative in deep contango"). It measures; it does **not** predict or opine (interpretation lives in the agent layer). This is the seam the rest of the Research OS hangs off.

## Dependencies
Depends only on `data/` (the HDB query layer `.data.hdb.*`) and `signals/` (seasonality, reused for the seasonal phase). It never imports `backtest/` or above — `.regime.breakdown` is a composition over plain inputs (PnL + labels), and `.strategy.commodityBT.__perf` is called lazily at runtime, not imported. Loads at import but never opens the HDB.

## Modules
- `regime.q` — `.regime.*`: the Market State Engine, the `regimes`-table build, and the regime-conditional breakdown.

## Key API
- `.regime.series[commodity;dates]` → a regime table (date + `curveState`/`volState`/`liqState`/`rollPhase`/`seasonPhase` labels + `slopePct`/`volPct`/`volumePct`).
- `.regime.label[commodity;date]` → the labels for one (commodity,date), computed causally over history up to that date.
- `.regime.breakdown[pnlByDate;regimeLabels;axis]` → per-regime-bucket performance (reusing `.strategy.commodityBT.__perf`) + a `(blended)` all-data row that reproduces the backtest headline.
- `.regime.buildTable[commodities;hdbPath]` → builds the splayed `regimes` table; `.regime.open[hdbPath]` loads it (via `get`, not `\l`, so it never changes the working directory).

## Axes
- `curveState` — front-vs-deferred relative slope: `backwardation` / `contango` / `flat`, with `slopePct` (rolling percentile vs its own history).
- `volState` — realized vol of front-settle log returns (roll-day returns zeroed), annualised → `low`/`normal`/`high` by `volPct`.
- `liqState` — front-volume rolling percentile → `thin`/`normal`/`deep`.
- `rollPhase` — front days-to-expiry bucket → `near`/`mid`.
- `seasonPhase` — `peak`/`trough`/`neutral` from `signals/seasonality`.

## Notes
- Thresholds live in `.cfg.regime` (no magic numbers). Percentiles are causal trailing-window fractions (no look-ahead, no nearest-rank index rounding).
- Pure core (`__panelFromLong` + `__labelPanel`) is HDB-free and known-answer tested on synthetic data (`tests/regime/`).
- The `regimes` table is splayed alongside `futures` under `.cfg.paths.hdb` (gitignored); build with `q scripts/build_regimes.q` after `ingest_hdb.q`. Demo: `apps/examples/regime_conditional_backtest.q`.
- Additive layer: the engine/backtest are untouched, the full suite stays byte-identical.
