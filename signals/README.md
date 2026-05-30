# signals/ — Alpha signal library (ARCHITECTURE.md §1, §6)

## Purpose
Deterministic signal construction for the strategies — currently the seasonality overlay.

## Dependencies
Depends on `core/` (and `models/` for curve types). Used by `backtest/` and `calibration/` (seasonal de-trending).

## Modules
- `seasonality.q` — `.commodity.seasonality.*`: the seasonal log-forward overlay (`factor[timeYears;seasonCfg]`, smooth annual cosine or explicit 12-vector) and `fitMonthlyFactors` (recovers a 12-vector seasonal pattern from a demeaned curve panel).

## Key API
`.commodity.seasonality.factor`, `.commodity.seasonality.fitMonthlyFactors`.

## Notes
- The signal-building half of the commodity harness still lives in `backtest/commodityStrategies.q` (`.strategy.path.commoditySignals`); splitting it here is a later step.
- Most relevant to gas/power; crude is only mildly seasonal.
