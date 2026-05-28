# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

qFDM is a kdb+/q pricing and risk framework. The equity FDM core (Black-Scholes, Crank-Nicolson, American, barriers, local vol) is the *technical validation case*; the long-term target is commodity (oil, power/electricity) and multi-asset options. Treat AAPL/equity tests as proof-of-correctness, not the end goal.

Current version: see `.qfdm.version` set at the bottom of `lib/init.q` (today: 0.41).

Test status: `q tests/run_all_tests.q` → **260 passed / 0 failed**.

v0.35.2 strengthened `.commodity.mrjump` validation with four benchmark tests (lambda-zero call/put vs Schwartz closed-form, stationary-variance limit, jump-intensity monotonicity, jump-mean sensitivity). No pricing logic was changed.

v0.36 adds `.commodity.modelreport` for cross-model commodity option comparison and scenario PnL across Black-76, Schwartz one-factor, Schwartz two-factor, and mean-reverting jump models. Reporting layer only — no pricing-formula changes.

v0.37 adds commodity model Greeks/sensitivities (finite-difference) inside `.commodity.modelreport` across all four commodity models, plus a union-joined `greeksAllModels` table, a `greeksSummary` aggregator, and a `runComparisonWithGreeks` convenience wrapper. Common random numbers (shared `mcConfig.randomSeed`) are reused across mrjump base/up/down to reduce FD noise. Reporting/sensitivity layer only — no pricing-formula changes.

v0.38 adds commodity model-disagreement risk metrics and alerting inside `.commodity.modelreport`: `primarySensitivity`, `priceDisagreement`, `greeksDisagreement`, `scenarioDisagreement`, `disagreementAlerts`, `modelDisagreementReport`, and the top-level `runComparisonRisk` wrapper. Status policy: zero OK rows -> `ERROR; 1..(minimumOkModels-1) -> `warning; otherwise `OK. Reporting-layer extension only — no pricing-formula changes.

v0.39 adds portfolio-level commodity model-disagreement aggregation inside `.commodity.modelreport`: `validatePortfolioPositions`, `runPositionRisk` (tagged per-position tables), `runPortfolioRisk` (uj-merged book of priceRows/greeksRows/scenarioPnlRows/alertRows + positionSummary), `portfolioAlertSummary`, `portfolioDisagreementExposure`, `worstOffenders` (top-N by metric — column is `offenderRank` to avoid shadowing q's `rank`), `portfolioRiskDashboard`, and `runPortfolioDashboard`. Per-position try-catch isolates failures so one bad position never crashes the book. Reporting/aggregation layer only — no pricing-formula changes.

v0.40 adds commodity model-risk limit monitoring on top of the portfolio disagreement outputs: `defaultModelRiskLimitConfig`, `checkLimit` (OK/warning/breach/ERROR), `checkPortfolioModelRiskLimits` (seven limits: grossPriceRangeExposure, maxScenarioPnlRange, maxPrimarySensitivityRange, maxVolatilityVegaRange, maxJumpSensitivity, warningAlertCount, errorTradeCount), `limitStatusSummary` (priority ERROR > breach > warning > OK), `limitBreachReport` (sorted, with `breachRank` column — `rank` would shadow q's keyword), and `runPortfolioRiskWithLimits` wrapper. All inputs are float-cast inside `checkLimit` so the resulting table has consistent column types regardless of whether a metric is a long count or a float exposure. Reporting/risk-control layer only — no pricing-formula changes.

v0.41 adds in-memory time-series limit history and trend analytics inside `.commodity.modelreport`: `limitSnapshot` (tags a limit table with `runId`/`runDate`/`portfolioName`), `limitSummarySnapshot`, `appendLimitHistory`, `limitHistorySummary` (per-limit observation/OK/warning/breach/error counts plus rates and latest), `limitBreachTrend` (latest vs avg-of-lookback, with worsening/improving/flat direction), `repeatedBreaches`, `limitHistoryDashboard`, and the top-level `runPortfolioRiskWithLimitHistory` wrapper. The wrapper bundles run metadata and existing histories into dicts because q lambdas allow at most 8 parameters. Reporting/aggregation layer only — no persistence and no pricing-formula changes.

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
