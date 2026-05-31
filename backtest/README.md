# backtest/ — Backtest layer (ARCHITECTURE.md §1)

## Purpose
The generic strategy engine + registry, the equity and commodity strategy suites, the signal-augmented-path harness, walk-forward, scoring, and cross-commodity robustness.

## Dependencies
Depends on `core/`, `models/`, `calibration/`, `signals/`, `data/`, and `execution/`. The top of the compute stack; called by `apps/`. The replay engine (`replay.q`) additionally composes `state/` (R9), `curve/` (R10), and `roll/` (R11) — all downward — and opens nothing at import.

## Modules
- `strategy.q` — `.strategy.*`: the generic registry-based engine (`register`/`run`/`runAndSummarize`/`defaultConfig` dispatcher), the path adapters (`path.fromFutures`/`fromCorrelated`/`fromFuturesCurve`/`ensemble`/…), shared hedge accounting (`__hedgeStep`/`__portfolioValue`), and the equity strategy suite (gammaScalp, shortVariance, calendarRoll, riskReversal, modelDisagreement, deltaVegaHedge, longVol, collarTailHedge, putRatioBackspread, ironCondor, barrierHedge, jumpPremium, dispersion, powerSpikeCapture, commodityCalendar, sparkSpread, crackSpread).
- `commodityStrategies.q` — `.strategy.commodityBT.*` + `.strategy.path.commoditySignals`: the commodity backtest core (`coreInit`/`coreStep`/`coreSummary` — wired through `.exec.fill`), the signal-augmented path builder, the 7 commodity strategies (convenienceYieldCarry, chiReversion, timeSeriesMomentum, twoTimescale, storageCashCarry, carryMomentumCombo, curveRelativeValue), `runSuite`, `walkForward`, `crossCommodity`.
- `replay.q` — `.backtest.replay.*` (**Research OS R12, the KEYSTONE**): the **event-driven replay engine** — a deterministic date-by-date fold that composes the as-of door (R9), the curve (R10), and roll discipline (R11) into realistic, cost-and-roll-aware PnL on the **actual contracts**, plus an auditable run record (R13 consumes it). A **NEW path** that REUSES the unchanged `.exec.fill` and research-mode engine (it does not modify them).

## Key API
`.strategy.register`, `.strategy.run`/`runAndSummarize`, `.strategy.defaultConfig`, `.strategy.path.commoditySignals`, `.strategy.commodityBT.runSuite`/`walkForward`/`crossCommodity`.
- `.backtest.replay.run[strategy;comm;fromDate;toDate;cfg]` → an auditable run record `` `meta`steps`rollEvents`provenance ``: per-step it builds the as-of Market State (R9), decides the active contract as of the prior date (R11's `days_before_expiry` engine, the causally-honest convention that tracks the research `heldSeq`), re-derives the as-of return on that contract, generates orders (target − current), fills via `.exec.fill` (+ opt-in realism), books PnL (price − cost − financing), and logs everything. Deterministic (no RNG).
- `.backtest.replay.faithfulness[strategy;comm;tol]` — the **look-ahead detector**: a frictionless replay reproduces the research-mode PnL to a tight tolerance on every date where the roll map names the same contract; a divergence there is look-ahead or a roll/fill bug. The few trailing **data-edge** dates (research holds a delisted contract) are reported as `nEdge`/`edgeDates`, not hidden.
- `.backtest.replay.persist`/`open` — thin opt-in write of the run's `steps` to the splayed, **unpartitioned** `replayRuns` (gitignored; tests synthetic).

## Notes
- Both engines book **net-of-execution** PnL via the `execution/` layer (`.exec.fill`); the default exec config is frictionless and byte-identical to the legacy flat cost. The commodity BT routes through `coreInit`/`coreStep` (v0.60); the equity engine routes its delta-hedge leg through `.strategy.__hedgeStep`/`__hedgeInit` (v0.61) and **every option-leg trade** (entry/roll/book/vega/lifecycle) through `.strategy.__legCost` (v0.63, step 4c) — Step 4 is fully complete, `deltaPV==stepPnl` preserved throughout.
- Per-strategy default configs live in `.cfg.strategy.*`; pass a per-run `exec` sub-config to enable slippage / cap realism (threaded into the hedge via `.strategy.__withExec` for all hedged strategies; option-leg slippage reaches `__legCost` via `stratCfg` directly).
- `commodityBT.coreSummary` reports GROSS vs NET (grossSharpe) + cost attribution; the walk-forward `aggregate`/`splits` schemas are fixed (tests pin them).
