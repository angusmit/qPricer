# calibration/ — Calibration layer (ARCHITECTURE.md §1)

## Purpose
Fit model parameters to market data: implied vol, vol surface, the generic optimizer/objective, single-curve Schwartz calibration, the Schwartz-Smith Kalman-filter MLE, and model-quality reports.

## Dependencies
Depends on `core/` and `models/`. Called by `analytics/` and `backtest/` (the commodity signal harness).

## Modules
- `iv.q` — `.iv.*`: implied-volatility solver + option-chain IVs.
- `surface.q` — `.surface.*`: volatility-surface construction + surface-based pricing.
- `objective.q` — `.objective.*`: objective functions (e.g. `rmse`) for fitting.
- `calibration.q` — `.calibration.*`: the generic calibration/search helpers (`bestCalibrationResult`).
- `calibrateCurve.q` — `.commodity.curveCal.*` + `.commodity.calibrateCurve`: fit `schwartz2` to a (tenor,price) curve; per-date convenience-yield series.
- `kalmanSchwartzSmith.q` — `.commodity.kalman.*`: the two-factor Schwartz-Smith state-space Kalman filter + MLE (`estimate`/`filter`/`panelFromCurveHistory`).
- `calibreport.q` / `modelcheck.q` / `modelcompare.q` — `.calibreport`/`.modelcheck`/`.modelcompare`: calibration reports + model-quality / consistency comparison.

## Key API
`.iv.impliedVolatility`, `.commodity.calibrateCurve[marketCurve;modelType;calCfg]`, `.commodity.curveCal.convenienceYieldSeries`, `.commodity.kalman.estimate`/`filter`.

## Notes
- Calibration bounds / grid-step defaults are externalized to `.cfg.calib.*` (curve / curveSeries / kalmanParams / kalmanEst); default config is byte-identical to the prior hardcoded values.
- Schwartz models are log-price: non-positive prices (e.g. the April-2020 WTI negatives) raise a controlled error.
