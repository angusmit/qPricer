# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

qFDM is a kdb+/q pricing and risk framework. The equity FDM core (Black-Scholes, Crank-Nicolson, American, barriers, local vol) is the *technical validation case*; the long-term target is commodity (oil, power/electricity) and multi-asset options. Treat AAPL/equity tests as proof-of-correctness, not the end goal.

Current version: see `.qfdm.version` set at the bottom of `lib/init.q` (today: 0.48).

Test status: `q tests/run_all_tests.q` → **322 passed / 0 failed**.

### Test runners (v0.48)
- `q tests/run_all_tests.q` — full suite + timing report (slowest-20, per-group subtotals).
- `q tests/run_group.q <group> -q` — single group, e.g. `strategy`, `core`, `commodity`.
- `q tests/run_fast_tests.q` — fast tier (deterministic core/greeks/market/iv/infra), < ~2 min.
- `q tests/run_smoke_tests.q` — very short core smoke, < 30 s.

All runners reuse `.test.suites` from `run_all_tests.q` via a `.test.skipAutoRun` flag — no duplication.

v0.35.2 strengthened `.commodity.mrjump` validation with four benchmark tests (lambda-zero call/put vs Schwartz closed-form, stationary-variance limit, jump-intensity monotonicity, jump-mean sensitivity). No pricing logic was changed.

v0.36 adds `.commodity.modelreport` for cross-model commodity option comparison and scenario PnL across Black-76, Schwartz one-factor, Schwartz two-factor, and mean-reverting jump models. Reporting layer only — no pricing-formula changes.

v0.37 adds commodity model Greeks/sensitivities (finite-difference) inside `.commodity.modelreport` across all four commodity models, plus a union-joined `greeksAllModels` table, a `greeksSummary` aggregator, and a `runComparisonWithGreeks` convenience wrapper. Common random numbers (shared `mcConfig.randomSeed`) are reused across mrjump base/up/down to reduce FD noise. Reporting/sensitivity layer only — no pricing-formula changes.

v0.38 adds commodity model-disagreement risk metrics and alerting inside `.commodity.modelreport`: `primarySensitivity`, `priceDisagreement`, `greeksDisagreement`, `scenarioDisagreement`, `disagreementAlerts`, `modelDisagreementReport`, and the top-level `runComparisonRisk` wrapper. Status policy: zero OK rows -> `ERROR; 1..(minimumOkModels-1) -> `warning; otherwise `OK. Reporting-layer extension only — no pricing-formula changes.

v0.39 adds portfolio-level commodity model-disagreement aggregation inside `.commodity.modelreport`: `validatePortfolioPositions`, `runPositionRisk` (tagged per-position tables), `runPortfolioRisk` (uj-merged book of priceRows/greeksRows/scenarioPnlRows/alertRows + positionSummary), `portfolioAlertSummary`, `portfolioDisagreementExposure`, `worstOffenders` (top-N by metric — column is `offenderRank` to avoid shadowing q's `rank`), `portfolioRiskDashboard`, and `runPortfolioDashboard`. Per-position try-catch isolates failures so one bad position never crashes the book. Reporting/aggregation layer only — no pricing-formula changes.

v0.40 adds commodity model-risk limit monitoring on top of the portfolio disagreement outputs: `defaultModelRiskLimitConfig`, `checkLimit` (OK/warning/breach/ERROR), `checkPortfolioModelRiskLimits` (seven limits: grossPriceRangeExposure, maxScenarioPnlRange, maxPrimarySensitivityRange, maxVolatilityVegaRange, maxJumpSensitivity, warningAlertCount, errorTradeCount), `limitStatusSummary` (priority ERROR > breach > warning > OK), `limitBreachReport` (sorted, with `breachRank` column — `rank` would shadow q's keyword), and `runPortfolioRiskWithLimits` wrapper. All inputs are float-cast inside `checkLimit` so the resulting table has consistent column types regardless of whether a metric is a long count or a float exposure. Reporting/risk-control layer only — no pricing-formula changes.

v0.41 adds in-memory time-series limit history and trend analytics inside `.commodity.modelreport`: `limitSnapshot` (tags a limit table with `runId`/`runDate`/`portfolioName`), `limitSummarySnapshot`, `appendLimitHistory`, `limitHistorySummary` (per-limit observation/OK/warning/breach/error counts plus rates and latest), `limitBreachTrend` (latest vs avg-of-lookback, with worsening/improving/flat direction), `repeatedBreaches`, `limitHistoryDashboard`, and the top-level `runPortfolioRiskWithLimitHistory` wrapper. The wrapper bundles run metadata and existing histories into dicts because q lambdas allow at most 8 parameters. Reporting/aggregation layer only — no persistence and no pricing-formula changes.

v0.42 adds a generic, registry-based strategy/backtest engine (`.strategy`) with data-source-agnostic path adapters (synthetic GBM + Barchart-normalised) and a first delta-hedged gamma-scalping strategy (`.strategy.gammaScalp`) that self-registers at load. The driver `.strategy.run` contains no strategy-specific branching; it folds `stepFn` over path rows via scan and builds the result table column-wise once. Gamma scalping produces per-step P&L attribution (option, hedge, financing, theta, theoretical-gamma) and a `gammaReconResidual` cross-check against the existing greeks engine — within 1% of initial option price on small-step paths. Calls `.engine.priceOption` and `.greeks.calculateGreeks`; no pricing-formula or routing changes.

v0.43 adds a short-variance (variance-risk-premium) strategy (`.strategy.shortVariance`) on the generic strategy engine, sharing a single delta-hedge accounting helper (`.strategy.__hedgeStep` + `.strategy.__hedgeInit`) with gamma scalping, and proves registry extensibility on a real second strategy. shortVariance builds a two-leg ATM straddle from the input trade, gates entry on `impliedVol > forecastVol + entryMargin`, collects premium up front, and runs delta-hedged with the same accounting machinery. The driver is untouched. Gamma scalp output is byte-identical post-refactor (totalPnl, intervalRebals, residual all match v0.42 PASS lines). No pricing-formula changes.

v0.44 adds a calendar-roll strategy (`.strategy.calendarRoll`) on the generic engine — the first with a mutable leg set and roll-event accounting via a portfolio-value identity — validating the registry and driver for path-dependent position lifecycle without engine changes. The leg book is a q TABLE held in state; expiring legs are settled at intrinsic from payoff (never priced via the FDM grid at zero expiry); replacement legs are opened at fresh tenors. Accounting uses `stepPnl == positionPnl + rollPnl + hedgePnl + financingPnl − txnCost` (max identity residual = 0 in the accounting test). gammaReconResidual is computed over non-roll steps only; rolls are reported separately. `.strategy.__hedgeStep` is reused for the hedge leg only; its `positionPnl` / `stepPnl` are discarded because the position mutates. Gamma scalp + short variance outputs unchanged. No pricing-formula changes.

v0.45 completes the strategy suite (`.strategy.riskReversal`, `.strategy.modelDisagreement`, `.strategy.deltaVegaHedge`) and adds a common-random-number multi-path ensemble portfolio runner with a strategy-native performance + correlation dashboard. Added `.strategy.__portfolioValue` helper (PV = cash + legMarkSum + hedgePosition*spot). Three new strategies self-register via the registry. New `.strategy.path.ensemble` for common random numbers across strategies. `.strategy.portfolio.runEnsemble`, `.performanceByStrategy` (mean/std/percentiles/winRate/sharpeLike), `.strategyCorrelation` (N×N matrix), `.dashboard`. Each strategy's accounting test uses `__portfolioValue` to compute deltaPV from state components independently and asserts equality to `stepPnl` from the result table (independent revaluation, not a tautology — max residuals ≤ 2e-14 across the three new strategies and the rewritten calendar test). Existing gamma scalp, short variance, calendar roll outputs unchanged. No pricing-formula changes.

v0.46 (Part A) adds three more equity-variant strategies on the generic engine: `.strategy.longVol` (mirror of shortVariance — buy ATM straddle when forecastVol > implied vol + entry margin), `.strategy.collarTailHedge` (collar mode: long underlying + long OTM put + short OTM call; tailHedge mode: long OTM puts sized to a premium-budget percent of notional, held outright with no external delta hedge), and `.strategy.putRatioBackspread` (short 1 near-ATM put + long ratioN OTM puts at a lower strike, optional delta hedge). All three self-register, reuse `.strategy.__hedgeStep` for the hedge leg where applicable, and ship accounting tests that compute deltaPV independently via `.strategy.__portfolioValue` and assert equality to result-table stepPnl within 1e-8 (max residuals ≤ 7e-15 across the three new strategies). Existing 9 strategies + ensemble dashboard byte-identical. Parts B–E (defined-risk + barrier; jump premium; correlated multi-asset dispersion; commodity futures-curve adapter + powerSpike + commodityCalendar) intentionally deferred to keep the milestone delivery whole rather than half-done. No pricing-formula changes; no engine/driver/registry/shared-helper changes.

v0.48 is a test-infrastructure milestone — no production code changes. Adds per-test timing to the runner (slowest-20 + per-group subtotals via `.z.p` deltas), a single-group runner (`tests/run_group.q <group>`) that reuses `.test.suites` via a `.test.skipAutoRun` flag, a fast-tier runner (`tests/run_fast_tests.q`), and a top-of-file doc-comment block listing the new commands. Tuned twelve non-canonical accounting/structural tests by coarsening the FDM grid (80×150 → 40×50) and shortening synthetic paths (7-9 → 5 steps); the deltaPV-vs-stepPnl identity holds at machine epsilon regardless of grid resolution, so tolerances were not relaxed. Canonical byte-identical reference tests (`gamma_scalp_synthetic` `totalPnl=0.6521844`, `short_variance_synthetic` `premium=11.97994`, all recon and gate tests) were NOT touched. Full suite wall time: **354 s** (~5.9 min, comfortably under the 10-min ceiling); strategy group 145 s, risk 95 s, stress 34 s, portfolio 28 s, commodity 15 s, all others < 12 s each. Slowest remaining are reference/recon tests that must keep their reference numbers. Trap recorded in `SKILL.md` (extends note 12): a bare `/` line — even one separating documentation comments — opens a multi-line comment block via `\l` and silently disables everything below.

v0.47 completes the v0.46 deferred work in five strategies + two new path adapters, lifting the suite from 9 to 14 strategies. Part B: `.strategy.ironCondor` (4 legs — short OTM put + long further-OTM put + short OTM call + long further-OTM call; held outright with optional delta hedge; analytic max-profit = netCredit, max-loss = widestSpreadWidth − netCredit, breakevens = shortStrike ± netCredit asserted to 1e-12), `.strategy.barrierHedge` (long upAndOut/downAndOut option from product trade fields, delta-hedged via `__hedgeStep`; on knock-out event the option is settled at rebate, hedge closed, position flat for remaining steps; tests use explicit FDM as Crank-Nicolson doesn't support barriers). Part C: `.strategy.jumpPremium` is the structural-model counterpart to modelDisagreement — at entry it prices the option under Merton (series, fast) or Bates (MC, tiny pathCount in tests) AND under BS; jumpPremium = jumpPx0 − bsPx0, gated by absolute/relative threshold; sells when market underprices jump risk (jumpPx0 > bsPx0), forceLong/forceShort overrides; marks on BS along the path for speed. Part D: new path adapter `.strategy.path.fromCorrelated` (Cholesky-correlated multi-asset GBM via `.correlation.__cholesky`, returns per-name path tables + weighted index path), and `.strategy.dispersion` (short index ATM straddle + long weighted constituent ATM straddles, per-underlying hedge book as a TABLE keyed by name, index leg held unhedged for correlation exposure; backs out implied avg correlation from `σ_idx² = w'Σw` and computes realized via pairwise log-return correlation; identity recovery asserted to 1e-10). Part E: new path adapter `.strategy.path.fromFuturesCurve` (`simple` or `mrjump` evolution of the front level + additive contango per tenor; returns frontPath standard-schema + long-format curveSnapshots + jumpCountsAtStep), `.strategy.powerSpikeCapture` (long Black-76 call on the front, gated by jumpCount > 0 OR deviation from mean > threshold·σ; on opening from flat the position-value jump is offset by cash outflow so stepPnl reports only `financing − txnCost`), `.strategy.commodityCalendar` (long near + short far futures on the curve; rolls near leg when elapsedYears ≥ nearTenor, opens next tenor; reuses calendarRoll's leg-table-in-state pattern). One additive multi-asset helper `.strategy.__portfolioValueMulti[cash;legMarkSum;hedgePositions;spots]` for dispersion's accounting test; the single-asset `.strategy.__portfolioValue` is byte-identical. Each new strategy ships an accounting test computing deltaPV from state via the appropriate portfolio-value helper and asserting equality to result-table stepPnl within 1e-8 (max residuals ≤ 7e-15 across all five new strategies). 9 existing strategies + ensemble dashboard byte-identical. No pricing-formula changes; no engine/driver/registry/shared-hedge-helper signature changes.

AAPL / Barchart data under `data/barchart/aapl/options_history/` is a **technical validation dataset only** — it exists to exercise the parser/backtest path on real-shaped CSVs, not because equities are the target asset class. The long-term target remains commodity (oil, power/electricity) and multi-asset options.

## Common commands

All commands assume you are in the repo root and invoke a q binary directly (this codebase has no build step — q loads `.q` files at runtime).

```
q tests/run_smoke_tests.q       / fast core sanity tests
q tests/run_all_tests.q         / full suite (17 grouped suites, hundreds of tests)
q tests/run_full_tests.q        / wrapper around run_all_tests.q
q benchmarks/run_benchmarks.q   / perf benchmarks, separate from unit tests
q tests/run_commodity_demo.q    / commodity stack demo (futures curve, Black-76, spreads, electricity)
q tests/run_barchart_backtest.q / Barchart historical option replay
q tests/run_daily_pricing.q     / daily risk orchestration
q tests/run_stress_test.q       / stress harness
q examples/<file>.q             / standalone scenarios (smoke_test_*, calculate_greeks, price_portfolio, ...)
```

### Running a single test

Tests are standalone scripts — each one starts with `\l lib/init.q` and can be run directly:

```
q tests/core/test_european_call.q
```

Inside the bulk runners (`run_all_tests.q`, `run_smoke_tests.q`, `run_benchmarks.q`) the harness reads each test file, **strips `\l ...` lines**, and `value`s the remaining source under a single already-loaded `lib/init.q`. Consequences when writing a new test:

- Top-level `\l lib/init.q` at the start of a test file is **stripped** by the runner — only direct `q <file>` invocations actually load init from inside the test.
- Anything other than `\l ...` runs in the runner's global scope; assume the lib is already loaded.
- A test passes if `value`ing the source does not throw; signal failure with `'"reason"`. There's no assertion framework — tests print PASS/FAIL and `'"..."` on failure.
- To add a new test to the suite, append the path to the relevant `.test.<group>Files` list in `tests/run_all_tests.q`.

## Architecture

### Loading pipeline

`lib/init.q` loads every module in dependency order and sets `.qfdm.loaded:1b` and `.qfdm.version`. When adding a new module, register it in `lib/init.q` — there is no auto-discovery.

### Pricing pipeline (equity FDM)

```
trade dict + marketData dict + model dict + config dict
    → .engine.__validateInputs
    → (if knock-in)  .engine.__priceKnockInViaParity         / knock-in = vanilla − knock-out
    → .engine.__runSolver        / dispatches to explicit / Crank-Nicolson, calls grid/payoff/boundary
    → .engine.__interpolatePriceFromGrid
    → returns dict: tradeId, underlying, optionType, unitPrice, notionalPrice, method, modelName
```

Layered modules (each one input-and-output, no hidden state):

- `product.q` — trade dict schema and validation; also `isKnockIn`, `knockInToKnockOut`.
- `market.q` / `marketbook.q` — single-symbol and multi-symbol market data.
- `model.q` — model selectors (Black-Scholes, local vol, Heston, Merton, Bates, SABR, Black-76, …).
- `grid.q`, `payoff.q`, `boundary.q` — solver building blocks.
- `solver.q` — explicit FDM and Crank-Nicolson backward induction.
- `engine.q` — **the public pricing API**: `.engine.priceOption`, `.engine.priceOptionWithGrid`.
- `greeks.q`, `risk.q`, `american.q`, `portfolio.q` — analytics built on top of `engine`.

### Beyond equity FDM

The library has grown well past the README's stated scope. Additional surfaces:

- **Monte Carlo stack** — `montecarlo.q`, `asian.q`, `correlation.q`, `basket.q`, `lookback.q`, `variance.q`, `pathdiag.q`.
- **Model calibration / comparison** — `iv.q`, `surface.q`, `calibration.q`, `objective.q`, `modelcheck.q`, `modelcompare.q`, `calibreport.q`, `convergence.q`.
- **Risk engine** — `riskdist.q`, `var.q`, `histscen.q`, `replay.q`, `limits.q`, `limitcheck.q`, `limitreport.q`, `dailyrisk.q`, `dashboard.q`.
- **Commodity / multi-asset** — `assetclass.q` (routing registry), `commodityFutures.q`, `commodityBlack76.q`, `schwartz.q` (`.commodity.schwartz`, one-factor log mean reversion, v0.33), `schwartz2.q` (`.commodity.schwartz2`, two-factor: spot + stochastic convenience yield, v0.34), `meanRevertingJump.q` (`.commodity.mrjump`, OU log-price + compound-Poisson jumps for spike/shock dynamics, v0.35), `commodityModelReport.q` (`.commodity.modelreport`, cross-model comparison + scenario PnL, v0.36), `commoditySpread.q`, `electricity.q`, `commodityBacktest.q`.
- **Real-data path** — `parser.q` (Barchart CSV loader — **deliberately standalone**, no qFDM deps) → `backtest.q`.
- **Reporting / infra** — `report.q`, `pnl.q`, `audit.q`, `regression.q`, `batch.q`, `result.q`, `timing.q`, `testutil.q`, `stress.q`, `perfdiag.q`, `perfopt.q`, `cache.q`.

### Asset class routing

`lib/assetclass.q` is the central registry mapping product type → asset class (`equity`, `commodity`, `electricity`, `rates`, `fx`) and model → model family. When adding a new product or model, update `.assetclass.__productMap` / `.assetclass.__modelMap` so routing stays correct.

## Conventions (enforced — `.claude/skills/q-pricing/SKILL.md` is the source of truth)

- **Namespaces are explicit and descriptive**: `.engine`, `.greeks`, `.commodity.futures`, `.parser.barchart`, etc. Never use `.d`, `.p`, `.x`, `.u`.
- **Private helpers use `.__`**: e.g. `.engine.__runSolver`, `.parser.barchart.__cleanCol`. Public functions sit at the top of the file; private helpers are below.
- **camelCase everywhere** — variables, functions, dict keys, table columns (`spotPrice`, `strikePrice`, `timeToExpiry`, `riskFreeRate`, `tradeId`).
- **Never shadow q built-ins** as names: not `value`, `count`, `select`, `key`, `var`, `type`, `string`, etc. Rename to project-specific terms (`optionValue`, `tradeCount`).
- **Short locals (`i`, `x`, `y`)** are fine only inside tiny anonymous expressions — never as public API, columns, or config keys.
- **Inputs are dicts** (trade, marketData, model, config). Outputs are dicts or tables with consistent column names. Don't introduce positional argument lists for public APIs.
- **Per-trade isolation in batch code** — portfolio/batch functions must wrap per-trade work so one bad trade returns `status=\`ERROR` rather than crashing the batch (see `portfolio.q`).
- **No silent pricing-logic changes.** A change that alters numerical output requires (a) a clear reason and (b) a test or recommended test, ideally with a benchmark comparison (closed-form, parity, limiting case).
- **Preserve working code** unless there is a clear bug or design problem. Prefer minimal patches over rewrites.

## q source cautions (silent failures)

These are the cautions enforced when writing or reviewing q in this repo. Each one parses or loads without an obvious error — the bug only surfaces later, which is why they are easy to miss.

- **No bare `/` separator lines.** A line consisting of just `/` opens a multiline comment block that runs until a standalone `\` line, silently commenting out everything in between (including function definitions). Use `/ text` on every comment line.
- **Avoid ambiguous `exp -1f`.** Tokenised as `exp - 1f` (subtraction), not negation. Always write `exp neg 1f` (or `exp neg expression`) when negating before `exp`/`log`/etc.
- **Typed numeric vector literals use a single trailing suffix.** Write `0.25 0.5 1 2 5f`, not `0.25 0.5 1f 2f 5f`. Patterns like `1 2 3f 4f 5f` or `1i 2i 3i` parse but do not produce the vector you expect.

## Things to know that aren't obvious from the file tree

- `parser.q` is intentionally standalone — do not add qFDM library calls into it; it is meant to be loadable independently for ingest jobs.
- `tests/run_all_tests.q` aborts via `'"Some tests failed: N failures"` if any test fails — don't catch this and continue.
- The `tests/realdata/test_barchart_*.q` suite is **fully synthetic** — those tests build option tables in q code and do not read any CSV. The CSV-backed runners (`tests/run_barchart_backtest.q`, `tests/realdata/run_barchart_backtest.q`) load every file under `data/barchart/aapl/options_history/` via `.parser.barchart.loadAll` and are **not part of `run_all_tests.q`** — they are demo/backtest runners and are sensitive to whatever CSVs sit in that directory.
- Knock-in barriers are priced as `vanilla − knock-out` (in-out parity) inside the engine — there is no dedicated knock-in solver.
- Local volatility is currently European-vanilla-only on the explicit FDM path; American / barrier + local-vol combinations are explicitly rejected.
- Crank-Nicolson currently supports European vanilla only.
- Portfolio multi-symbol pricing exists (`portfolio.q`, see multi-symbol tests) but assumes the caller supplies per-symbol market data.

## When to use the q-pricing skill

The repository ships `.claude/skills/q-pricing/SKILL.md`, which encodes the user's preferences for naming, namespacing, FDM/commodity/jump-diffusion/exotic checklists, code-review response format, and the long-term commodity roadmap. Reach for it when reviewing q code, designing a new pricer, evaluating a new technology, or answering maths-knowledge questions — those are the cases where its checklists carry real content beyond what's in this file.
