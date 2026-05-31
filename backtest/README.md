# backtest/ — Backtest layer (ARCHITECTURE.md §1)

## Purpose
The generic strategy engine + registry, the equity and commodity strategy suites, the signal-augmented-path harness, walk-forward, scoring, and cross-commodity robustness.

## Dependencies
Depends on `core/`, `models/`, `calibration/`, `signals/`, `data/`, and `execution/`. The top of the compute stack; called by `apps/`.

## Modules
- `strategy.q` — `.strategy.*`: the generic registry-based engine (`register`/`run`/`runAndSummarize`/`defaultConfig` dispatcher), the path adapters (`path.fromFutures`/`fromCorrelated`/`fromFuturesCurve`/`ensemble`/…), shared hedge accounting (`__hedgeStep`/`__portfolioValue`), and the equity strategy suite (gammaScalp, shortVariance, calendarRoll, riskReversal, modelDisagreement, deltaVegaHedge, longVol, collarTailHedge, putRatioBackspread, ironCondor, barrierHedge, jumpPremium, dispersion, powerSpikeCapture, commodityCalendar, sparkSpread, crackSpread).
- `commodityStrategies.q` — `.strategy.commodityBT.*` + `.strategy.path.commoditySignals`: the commodity backtest core (`coreInit`/`coreStep`/`coreSummary` — wired through `.exec.fill`), the signal-augmented path builder, the 7 commodity strategies (convenienceYieldCarry, chiReversion, timeSeriesMomentum, twoTimescale, storageCashCarry, carryMomentumCombo, curveRelativeValue), `runSuite`, `walkForward`, `crossCommodity`.

## Key API
`.strategy.register`, `.strategy.run`/`runAndSummarize`, `.strategy.defaultConfig`, `.strategy.path.commoditySignals`, `.strategy.commodityBT.runSuite`/`walkForward`/`crossCommodity`.

## Notes
- Both engines book **net-of-execution** PnL via the `execution/` layer (`.exec.fill`); the default exec config is frictionless and byte-identical to the legacy flat cost. The commodity BT routes through `coreInit`/`coreStep` (v0.60); the equity engine routes its delta-hedge leg through `.strategy.__hedgeStep`/`__hedgeInit` (v0.61) and **every option-leg trade** (entry/roll/book/vega/lifecycle) through `.strategy.__legCost` (v0.63, step 4c) — Step 4 is fully complete, `deltaPV==stepPnl` preserved throughout.
- Per-strategy default configs live in `.cfg.strategy.*`; pass a per-run `exec` sub-config to enable slippage / cap realism (threaded into the hedge via `.strategy.__withExec` for all hedged strategies; option-leg slippage reaches `__legCost` via `stratCfg` directly).
- `commodityBT.coreSummary` reports GROSS vs NET (grossSharpe) + cost attribution; the walk-forward `aggregate`/`splits` schemas are fixed (tests pin them).
