# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

qFDM is a kdb+/q pricing and risk framework. The equity FDM core (Black-Scholes, Crank-Nicolson, American, barriers, local vol) is the *technical validation case*; the long-term target is commodity (oil, power/electricity) and multi-asset options. Treat AAPL/equity tests as proof-of-correctness, not the end goal.

Current version: see `.qfdm.version` set at the bottom of `lib/init.q` (today: 0.32).

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
- **Commodity / multi-asset** — `assetclass.q` (routing registry), `commodityFutures.q`, `commodityBlack76.q`, `commoditySpread.q`, `electricity.q`, `commodityBacktest.q`.
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

## Things to know that aren't obvious from the file tree

- `parser.q` is intentionally standalone — do not add qFDM library calls into it; it is meant to be loadable independently for ingest jobs.
- `tests/run_all_tests.q` aborts via `'"Some tests failed: N failures"` if any test fails — don't catch this and continue.
- The Barchart fixtures under `data/barchart/aapl/options_history/` are required by `tests/realdata/test_barchart_*.q`; don't delete them when cleaning.
- Knock-in barriers are priced as `vanilla − knock-out` (in-out parity) inside the engine — there is no dedicated knock-in solver.
- Local volatility is currently European-vanilla-only on the explicit FDM path; American / barrier + local-vol combinations are explicitly rejected.
- Crank-Nicolson currently supports European vanilla only.
- Portfolio multi-symbol pricing exists (`portfolio.q`, see multi-symbol tests) but assumes the caller supplies per-symbol market data.

## When to use the q-pricing skill

The repository ships `.claude/skills/q-pricing/SKILL.md`, which encodes the user's preferences for naming, namespacing, FDM/commodity/jump-diffusion/exotic checklists, code-review response format, and the long-term commodity roadmap. Reach for it when reviewing q code, designing a new pricer, evaluating a new technology, or answering maths-knowledge questions — those are the cases where its checklists carry real content beyond what's in this file.
