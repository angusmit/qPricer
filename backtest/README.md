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
- The commodity backtest books **net-of-execution** PnL via the `execution/` layer (`.exec.fill`); the default exec config is frictionless and byte-identical to the legacy flat cost. The equity engine is not yet wired to execution (step 4b).
- Per-strategy default configs live in `.cfg.strategy.*`; pass a per-run `exec` sub-config to enable slippage / cap realism.
- `coreSummary` reports GROSS vs NET (grossSharpe) + cost attribution; the walk-forward `aggregate`/`splits` schemas are fixed (tests pin them).
