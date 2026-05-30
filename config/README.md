# config/ — Configuration layer (ARCHITECTURE.md §1, §2)

## Purpose
Holds the `.cfg` value files — every tunable library default as native q dicts, with optional per-environment overrides.

## Dependencies
None (pure value files). Loaded by `core/cfg.q` first, so `.cfg` is available to every layer.

## Modules
- `base.q` — **all** baseline defaults as domain-keyed dicts: `.cfg.fdm` (FDM grid), `.cfg.iv` (IV solver), `.cfg.mc` (Monte Carlo), `.cfg.calib.*` (curve/Kalman calibration), `.cfg.analytics.*` (model-report thresholds), `.cfg.paths` (data/output/HDB dirs), `.cfg.strategy.*` (every per-strategy config + `commoditySignals`), `.cfg.exec` (execution fill-and-cost). Each value equals the literal it replaced, with types and key-order preserved.
- `dev.q` — demo per-environment override (coarser FDM grid, smaller MC paths, an output path, a gammaScalp `deltaBand`). Loaded only when `QPRICER_ENV=dev`.

## Key API
`.cfg.*` (populated at startup). Override mechanism: `QPRICER_ENV=<env> q ...` (PowerShell: `$env:QPRICER_ENV='dev'; q ...`).

## Notes
- Unset `QPRICER_ENV` ⇒ base only ⇒ default/test behavior, byte-identical to the prior hardcoded values. An unknown env falls back to base with a warning.
- The config sweep (migration Step 2) is **complete** across core, calibration, analytics, paths, all backtest strategy/signal configs, and execution.
- Tabular reference data (calendars, instrument specs) belongs here too as CSV/keyed tables — none yet.
