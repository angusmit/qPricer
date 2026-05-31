# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

qFDM is a kdb+/q pricing and risk framework. The equity FDM core (Black-Scholes, Crank-Nicolson, American, barriers, local vol) is the *technical validation case*; the long-term target is commodity (oil, power/electricity) and multi-asset options. Treat AAPL/equity tests as proof-of-correctness, not the end goal.

Current version: see `.qfdm.version` set at the bottom of `core/init.q` (today: 0.63).

Test status: `q tests/run_all_tests.q` → **373 passed / 0 failed**.

**Layer layout (ARCHITECTURE.md §1).** Code is organized into layers, dependencies flowing downward only: `core/` (math/RNG/stats/infra + the loader `core/init.q` + the config loader `core/cfg.q`), `config/` (the `.cfg` value files: `base.q` + optional `{env}.q` overrides), `data/` (Barchart parser + the splayed HDB `.data.hdb`), `models/` (BS/FDM core + all pricers + pricing domain), `calibration/` (iv, surface, objective, calibrate-curve, Kalman MLE, model quality), `analytics/` (risk/VaR/scenarios/limits/portfolio/reporting/perf), `signals/` (seasonality), `execution/` (daily fill-and-cost simulation `.exec`, wired into the commodity backtest), `backtest/` (strategy engine + commodity suite + walk-forward), `scripts/` (HDB ingestion `ingest_hdb.q` + reserved CI/pipeline), `apps/` (examples + demos), `tests/` (flat). `portfolio/` (cross-strategy allocator) and `services/` (optional IPC) are reserved skeletons. Each dir has a `README.md` describing its modules + API. Load with `\l core/init.q`; run the suite with `q tests/run_all_tests.q`.

**Configuration (`.cfg`, v0.57–0.58, ARCHITECTURE.md §2 — Step 2 COMPLETE).** Tunable library defaults live in `.cfg` (domain-keyed dicts), populated at startup by `core/cfg.q`: it loads `config/base.q` (every value = the prior hardcoded literal) then, if `QPRICER_ENV` is set, `config/{env}.q` to override a subset via upsert. Unset `QPRICER_ENV` (the default/test path) loads base only, so behavior is byte-identical. Domains: `.cfg.fdm` (FDM grid/pricing), `.cfg.iv` (implied-vol solver), `.cfg.mc` (Monte Carlo defaults), `.cfg.calib.*` (curve + Kalman calibration bounds/grids/params), `.cfg.analytics.*` (model-report greek-bump / disagreement / model-risk-limit thresholds), `.cfg.paths` (default data/output dirs), and `.cfg.strategy.*` (v0.58 — every per-strategy `defaultConfig`, 17 equity + 7 commodity, plus `.cfg.strategy.commoditySignals` = the commodity signal-path config). Each module's default-config function now just returns its `.cfg.*` dict. To run with overrides: `QPRICER_ENV=dev q tests/...` (PowerShell: `$env:QPRICER_ENV='dev'; q ...`); `config/dev.q` is a demo override (FDM grid, MC paths, a path, and a gammaScalp `deltaBand`). The config sweep is now complete across every layer; `models/` and `signals/` needed no sweep (params are caller-supplied). One inline numeric grid remains by design inside `.strategy.realCurveCalendarRoll` (a fixed Crank-Nicolson 40×50 driver grid, not a `defaultConfig` — left inline like a fixture).

**Data layer / HDB (`.data.hdb`, v0.59, ARCHITECTURE.md §3 — Step 3).** Real Barchart data moves off per-run CSV re-parsing onto a **splayed, unpartitioned** kdb+ HDB (`data/hdb`, from `.cfg.paths.hdb`; gitignored). Table `futures` (columns `commodity` sym/`p#`, `contractYM` long, `expiry`/`firstDate`/`date` dates, OHLC+settle floats, `volume` long) holds every commodity's cross-section. **Splayed not date-partitioned** deliberately: ~48k rows total, so ~1900 daily partition dirs would be slow + Windows-hostile; a single splay is far faster. Ingestion (`scripts/ingest_hdb.q` → `.data.hdb.ingest`) **reuses `.parser.futures.loadAll` verbatim** (so expiry=MAX/firstDate=MIN dates are derived exactly as the parser), enumerates syms with `.Q.en`, writes column files with `set`. Query layer `.data.hdb.curveAt` / `curveHistory` / `dates` materialise the parser-schema long table and feed it to the existing `.parser.crude.curveAt` / `curveHistory`, so output is **byte-identical to the parser by construction**. `.data.hdb.open` uses `get` (NOT `\l`) so it never changes the process working directory. The parser is UNCHANGED and is now the CSV-read stage feeding ingestion. `core/init.q` loads the query layer but never opens the HDB at import (library-load independence; a fresh checkout has no HDB). Build it with `q scripts/ingest_hdb.q`; rebuild is idempotent. Measured: CRUDE data-load **1884 ms (CSV parse) → 396 ms (HDB) ≈ 4.8× faster**; full build 49,078 rows / 1,866 dates in ~0.4 s. Real-data examples prefer the HDB when present and fall back to the parser otherwise. The synthetic test `tests/realdata/test_hdb_ingest_query_equivalence.q` proves ingest↔query equivalence with no real data.

**Execution layer (`.exec`, v0.60–0.63, ARCHITECTURE.md §5 — Step 4 FULLY COMPLETE).** `execution/execution.q` is a generic daily fill-and-cost simulation that turns a strategy's target position change (an order) into a simulated FILL (optional slippage + a liquidity participation cap) and explicit COST components, so the backtest books REALIZED (net) PnL and can report GROSS vs NET. `.exec.fill[order;ctx;cfg]` (pure, stateless) returns `filledQty`/`fillPrice`/`proportionalCost`/`slippageCost`/`fixedCost`/`totalCost`: full fill unless `participationCap` binds (`filled=min(|order|,cap*barVolume)`), proportional cost = `rate*|filled|`, slippage cost = `|filled|*refPrice*bps/1e4` (adverse), optional fixed + size/ADV impact. Config `.cfg.exec` is **frictionless by default**; `.exec.__resolve` injects the strategy's own txnCostRate as proportionalRate, so the **default reproduces the legacy flat-cost path byte-identically**. Both backtests are wired: the **commodity** BT (v0.60: `.strategy.commodityBT.coreInit`/`coreStep`; order = vol-targeted Δposition, refPrice = front price) and the **equity** engine (v0.61: the shared hedge helper `.strategy.__hedgeStep`/`__hedgeInit` — covering underlying-rebalancing turnover for every hedged strategy in one place; the hedge order is the dollar NOTIONAL with refPrice=1 and a fill-fraction, so the generic proportional cost reproduces the price-scaled legacy hedge cost `|hedgeTrade|*spot*rate` exactly and the `deltaPV==stepPnl` identity is preserved). A per-run `exec` sub-config is threaded into the hedge inputs via `.strategy.__withExec` (wired for gammaScalp + shortVariance; absent → frictionless → byte-identical for all). Under a cap the held hedge lags target (residual exposure). **v0.63 (step 4c) closed the remainder so EVERY equity trade routes through `.exec.fill`:** (1) the per-run `exec` sub-config is threaded into the hedge inputs of ALL hedged strategies (init+step `__withExec`); (2) every discrete OPTION-LEG cost (entry/roll/book-entry/vega-trade/lifecycle-settle, plus the spread/futures-leg and inline-hedge costs) routes through the shared `.strategy.__legCost[notional;stratCfg]` helper — order = the dollar notional, refPrice=1, so proportional cost reproduces the legacy per-leg cost EXACTLY and slippage = `|notional|*bps` is the **premium-scaled option bid-ask** (typically the dominant cost). Suite **373/0** (all prior tests byte-identical: gamma_scalp 0.6521844, short_variance 11.97994, all 17 accounting-identity residuals unchanged to the digit; + 3 equity execution tests). Both desks now fully wired (hedge + option legs); the participation cap applies to the continuous hedge, discrete legs fill fully. Gross-vs-net demos: `apps/examples/execution_gross_vs_net.q` (commodity, HDB), `execution_gross_vs_net_equity.q` (gamma-scalp hedge slippage + iron-condor option-leg slippage, synthetic).

**Portfolio allocator (`.alloc`, v0.62, ARCHITECTURE.md §6/§9 — Step 5).** `portfolio/portfolio.q` allocates across strategies from their net-of-execution return panels (N strategies × T dates). `.alloc.weights[returns;method;cfg]` supports `equalWeight` (1/N baseline), `inverseVol`, `minVariance` (Σ⁻¹1), `riskParity` (ERC via cyclical coordinate descent — **the DEFAULT**, covariance-only so it doesn't overfit), `maxSharpe` (Σ⁻¹μ), `meanVariance` ((1/λ)Σ⁻¹μ); constraints (`.cfg.alloc`): `longOnly`+`fullyInvested` (capped-simplex bisection projection), `weightCap`, `turnoverPenalty` (proximal shrink-to-prev), optional covariance `shrinkage`. `.alloc.backtest`/`.alloc.compare` do a CAUSAL walk-forward (weights from each split's TRAIN columns only, applied OOS; reuses `.strategy.commodityBT.__splits`/`__perf`) and rank methods by OOS Sharpe with equalWeight as the baseline — the honest "does optimization beat 1/N OOS?" table. ADDITIVE: consumes existing returns, touches no upstream code; the existing **370 tests stay byte-identical** (+ 2 new alloc tests = **372/0**). Distinct namespace: NOT analytics' `.portfolio.*` (option book) nor the ensemble `.strategy.portfolio.*` (cross-path dashboard). Demo: `apps/examples/portfolio_allocation.q` (real commodity net returns via the HDB).

### Test runners (v0.48)
- `q tests/run_all_tests.q` — full suite + timing report (slowest-20, per-group subtotals).
- `q tests/run_group.q <group> -q` — single group, e.g. `strategy`, `core`, `commodity`.
- `q tests/run_fast_tests.q` — fast tier (deterministic core/greeks/market/iv/infra), < ~2 min.
- `q tests/run_smoke_tests.q` — very short core smoke, < 30 s.

All runners reuse `.test.suites` from `run_all_tests.q` via a `.test.skipAutoRun` flag — no duplication.

Full version history: **CHANGELOG.md**.

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
q apps/examples/<file>.q        / standalone scenarios (smoke_test_*, calculate_greeks, price_portfolio, ...)
```

### Running a single test

Tests are standalone scripts — each one starts with `\l core/init.q` and can be run directly:

```
q tests/core/test_european_call.q
```

Inside the bulk runners (`run_all_tests.q`, `run_smoke_tests.q`, `run_benchmarks.q`) the harness reads each test file, **strips `\l ...` lines**, and `value`s the remaining source under a single already-loaded `core/init.q`. Consequences when writing a new test:

- Top-level `\l core/init.q` at the start of a test file is **stripped** by the runner — only direct `q <file>` invocations actually load init from inside the test.
- Anything other than `\l ...` runs in the runner's global scope; assume the lib is already loaded.
- A test passes if `value`ing the source does not throw; signal failure with `'"reason"`. There's no assertion framework — tests print PASS/FAIL and `'"..."` on failure.
- To add a new test to the suite, append the path to the relevant `.test.<group>Files` list in `tests/run_all_tests.q`.

## Architecture

### Loading pipeline

`core/init.q` loads every module in layer dependency order and sets `.qfdm.loaded:1b` and `.qfdm.version`. When adding a new module, register it in `core/init.q` — there is no auto-discovery.

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

`models/assetclass.q` is the central registry mapping product type → asset class (`equity`, `commodity`, `electricity`, `rates`, `fx`) and model → model family. When adding a new product or model, update `.assetclass.__productMap` / `.assetclass.__modelMap` so routing stays correct.

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

- **No bare `/` separator lines.** A line consisting of just `/` opens a multiline comment block that runs until a standalone `\` line, silently commenting out everything in between (including function definitions). Use `/ text` (e.g. `/ ----`) on every comment line. *Workflow habit (recurring bite): the blank divider line inside a new file's header doc-block is the usual culprit and loads without error — after writing/editing any `.q` file or example, grep `^\s*/\s*$` before the first run. Symptom: a just-defined `.ns.fn` errors with the function name as the signal (undefined), or a script produces no output.*
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

## Milestone working conventions (apply to every task)

Shell / process hygiene
- Run every q script as `q file.q </dev/null` (stdin EOF guarantees q exits even if a script forgets
  `exit 0;`); also end scripts with `exit 0;`.
- Quick checks run FOREGROUND. Only the full suite and the cross-commodity example may be backgrounded
  — capture to a file under `scratch/` and WAIT; never start the next command while one runs, never
  let a returning shell kill it.
- Verify no orphan q processes at the end (`Get-Process q`).

Cleanup / filesystem
- Do NOT run rm / rmdir / Remove-Item / any cleanup script (the deny list blocks them and workarounds
  hang). Leave all scratch under the gitignored `scratch/` for me to wipe.
- Do NOT use `mkdir` — kdb+ creates directories on write; a Windows `mkdir -p` spawns a stray `-p` folder.

Byte-identical discipline
- Canonical numbers stay EXACT unless a change is explicitly intended: gamma_scalp totalPnl=0.6521844,
  short_variance premium=11.97994, crude 63.27->57.68, cross-commodity momentum mean +0.7403373 / 0.75
  of cells, gas carry -0.0243 -> +0.1773. If a pinned number moves unexpectedly, it's a bug — fix it.
- Tests use SYNTHETIC data only; real CSVs + the HDB are gitignored, touched only by apps/examples + scripts.

Testing workflow
- During development run only the affected group (`q tests/run_group.q <group> -q </dev/null`). Run the
  full suite (`q tests/run_all_tests.q </dev/null`) only as the pre-commit gate.

Git / docs
- Never stage/commit without explicit approval; never stage bridge.q or test.py.
- Every milestone updates EVERY affected README (the touched layers' + the root README.md) so each
  reflects its dir's current modules/API — no partial/stale READMEs.

(q source-language rules — bare "/" lines, two-letter builtins, type fidelity — live in .claude/skills/q-pricing/SKILL.md.)