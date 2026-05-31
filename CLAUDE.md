# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

qFDM is a kdb+/q pricing and risk framework. The equity FDM core (Black-Scholes, Crank-Nicolson, American, barriers, local vol) is the *technical validation case*; the long-term target is commodity (oil, power/electricity) and multi-asset options. Treat AAPL/equity tests as proof-of-correctness, not the end goal.

Current version: see `.qfdm.version` set at the bottom of `core/init.q` (today: 0.69).

Test status: `q tests/run_all_tests.q` → **389 passed / 0 failed**.

**Layer layout (ARCHITECTURE.md §1).** Code is organized into layers, dependencies flowing downward only: `core/` (math/RNG/stats/infra + the loader `core/init.q` + the config loader `core/cfg.q`), `config/` (the `.cfg` value files: `base.q` + optional `{env}.q` overrides), `data/` (Barchart parser + the splayed HDB `.data.hdb`), `models/` (BS/FDM core + all pricers + pricing domain), `calibration/` (iv, surface, objective, calibrate-curve, Kalman MLE, model quality), `analytics/` (risk/VaR/scenarios/limits/portfolio/reporting/perf), `signals/` (seasonality), `execution/` (daily fill-and-cost simulation `.exec`, wired into the commodity backtest), `backtest/` (strategy engine + commodity suite + walk-forward), `scripts/` (HDB ingestion `ingest_hdb.q` + reserved CI/pipeline), `apps/` (examples + demos), `tests/` (flat). `portfolio/` (cross-strategy allocator) and `services/` (optional IPC) are reserved skeletons. Each dir has a `README.md` describing its modules + API. Load with `\l core/init.q`; run the suite with `q tests/run_all_tests.q`.

**Configuration (`.cfg`, v0.57–0.58, ARCHITECTURE.md §2 — Step 2 COMPLETE).** Tunable library defaults live in `.cfg` (domain-keyed dicts), populated at startup by `core/cfg.q`: it loads `config/base.q` (every value = the prior hardcoded literal) then, if `QPRICER_ENV` is set, `config/{env}.q` to override a subset via upsert. Unset `QPRICER_ENV` (the default/test path) loads base only, so behavior is byte-identical. Domains: `.cfg.fdm` (FDM grid/pricing), `.cfg.iv` (implied-vol solver), `.cfg.mc` (Monte Carlo defaults), `.cfg.calib.*` (curve + Kalman calibration bounds/grids/params), `.cfg.analytics.*` (model-report greek-bump / disagreement / model-risk-limit thresholds), `.cfg.paths` (default data/output dirs), and `.cfg.strategy.*` (v0.58 — every per-strategy `defaultConfig`, 17 equity + 7 commodity, plus `.cfg.strategy.commoditySignals` = the commodity signal-path config). Each module's default-config function now just returns its `.cfg.*` dict. To run with overrides: `QPRICER_ENV=dev q tests/...` (PowerShell: `$env:QPRICER_ENV='dev'; q ...`); `config/dev.q` is a demo override (FDM grid, MC paths, a path, and a gammaScalp `deltaBand`). The config sweep is now complete across every layer; `models/` and `signals/` needed no sweep (params are caller-supplied). One inline numeric grid remains by design inside `.strategy.realCurveCalendarRoll` (a fixed Crank-Nicolson 40×50 driver grid, not a `defaultConfig` — left inline like a fixture).

**Data layer / HDB (`.data.hdb`, v0.59, ARCHITECTURE.md §3 — Step 3).** Real Barchart data moves off per-run CSV re-parsing onto a **splayed, unpartitioned** kdb+ HDB (`data/hdb`, from `.cfg.paths.hdb`; gitignored). Table `futures` (columns `commodity` sym/`p#`, `contractYM` long, `expiry`/`firstDate`/`date` dates, OHLC+settle floats, `volume` long) holds every commodity's cross-section. **Splayed not date-partitioned** deliberately: ~48k rows total, so ~1900 daily partition dirs would be slow + Windows-hostile; a single splay is far faster. Ingestion (`scripts/ingest_hdb.q` → `.data.hdb.ingest`) **reuses `.parser.futures.loadAll` verbatim** (so expiry=MAX/firstDate=MIN dates are derived exactly as the parser), enumerates syms with `.Q.en`, writes column files with `set`. Query layer `.data.hdb.curveAt` / `curveHistory` / `dates` materialise the parser-schema long table and feed it to the existing `.parser.crude.curveAt` / `curveHistory`, so output is **byte-identical to the parser by construction**. `.data.hdb.open` uses `get` (NOT `\l`) so it never changes the process working directory. The parser is UNCHANGED and is now the CSV-read stage feeding ingestion. `core/init.q` loads the query layer but never opens the HDB at import (library-load independence; a fresh checkout has no HDB). Build it with `q scripts/ingest_hdb.q`; rebuild is idempotent. Measured: CRUDE data-load **1884 ms (CSV parse) → 396 ms (HDB) ≈ 4.8× faster**; full build 49,078 rows / 1,866 dates in ~0.4 s. Real-data examples prefer the HDB when present and fall back to the parser otherwise. The synthetic test `tests/realdata/test_hdb_ingest_query_equivalence.q` proves ingest↔query equivalence with no real data.

**Execution layer (`.exec`, v0.60–0.63, ARCHITECTURE.md §5 — Step 4 FULLY COMPLETE).** `execution/execution.q` is a generic daily fill-and-cost simulation that turns a strategy's target position change (an order) into a simulated FILL (optional slippage + a liquidity participation cap) and explicit COST components, so the backtest books REALIZED (net) PnL and can report GROSS vs NET. `.exec.fill[order;ctx;cfg]` (pure, stateless) returns `filledQty`/`fillPrice`/`proportionalCost`/`slippageCost`/`fixedCost`/`totalCost`: full fill unless `participationCap` binds (`filled=min(|order|,cap*barVolume)`), proportional cost = `rate*|filled|`, slippage cost = `|filled|*refPrice*bps/1e4` (adverse), optional fixed + size/ADV impact. Config `.cfg.exec` is **frictionless by default**; `.exec.__resolve` injects the strategy's own txnCostRate as proportionalRate, so the **default reproduces the legacy flat-cost path byte-identically**. Both backtests are wired: the **commodity** BT (v0.60: `.strategy.commodityBT.coreInit`/`coreStep`; order = vol-targeted Δposition, refPrice = front price) and the **equity** engine (v0.61: the shared hedge helper `.strategy.__hedgeStep`/`__hedgeInit` — covering underlying-rebalancing turnover for every hedged strategy in one place; the hedge order is the dollar NOTIONAL with refPrice=1 and a fill-fraction, so the generic proportional cost reproduces the price-scaled legacy hedge cost `|hedgeTrade|*spot*rate` exactly and the `deltaPV==stepPnl` identity is preserved). A per-run `exec` sub-config is threaded into the hedge inputs via `.strategy.__withExec` (wired for gammaScalp + shortVariance; absent → frictionless → byte-identical for all). Under a cap the held hedge lags target (residual exposure). **v0.63 (step 4c) closed the remainder so EVERY equity trade routes through `.exec.fill`:** (1) the per-run `exec` sub-config is threaded into the hedge inputs of ALL hedged strategies (init+step `__withExec`); (2) every discrete OPTION-LEG cost (entry/roll/book-entry/vega-trade/lifecycle-settle, plus the spread/futures-leg and inline-hedge costs) routes through the shared `.strategy.__legCost[notional;stratCfg]` helper — order = the dollar notional, refPrice=1, so proportional cost reproduces the legacy per-leg cost EXACTLY and slippage = `|notional|*bps` is the **premium-scaled option bid-ask** (typically the dominant cost). Suite **373/0** (all prior tests byte-identical: gamma_scalp 0.6521844, short_variance 11.97994, all 17 accounting-identity residuals unchanged to the digit; + 3 equity execution tests). Both desks now fully wired (hedge + option legs); the participation cap applies to the continuous hedge, discrete legs fill fully. Gross-vs-net demos: `apps/examples/execution_gross_vs_net.q` (commodity, HDB), `execution_gross_vs_net_equity.q` (gamma-scalp hedge slippage + iron-condor option-leg slippage, synthetic).

**Portfolio allocator (`.alloc`, v0.62, ARCHITECTURE.md §6/§9 — Step 5).** `portfolio/portfolio.q` allocates across strategies from their net-of-execution return panels (N strategies × T dates). `.alloc.weights[returns;method;cfg]` supports `equalWeight` (1/N baseline), `inverseVol`, `minVariance` (Σ⁻¹1), `riskParity` (ERC via cyclical coordinate descent — **the DEFAULT**, covariance-only so it doesn't overfit), `maxSharpe` (Σ⁻¹μ), `meanVariance` ((1/λ)Σ⁻¹μ); constraints (`.cfg.alloc`): `longOnly`+`fullyInvested` (capped-simplex bisection projection), `weightCap`, `turnoverPenalty` (proximal shrink-to-prev), optional covariance `shrinkage`. `.alloc.backtest`/`.alloc.compare` do a CAUSAL walk-forward (weights from each split's TRAIN columns only, applied OOS; reuses `.strategy.commodityBT.__splits`/`__perf`) and rank methods by OOS Sharpe with equalWeight as the baseline — the honest "does optimization beat 1/N OOS?" table. ADDITIVE: consumes existing returns, touches no upstream code; the existing **370 tests stay byte-identical** (+ 2 new alloc tests = **372/0**). Distinct namespace: NOT analytics' `.portfolio.*` (option book) nor the ensemble `.strategy.portfolio.*` (cross-path dashboard). Demo: `apps/examples/portfolio_allocation.q` (real commodity net returns via the HDB).

**Registry spine (`.registry` / `.contracts`, v0.68, ARCHITECTURE.md §13-R2).** The plug-in spine: `core/registry.q` is a GENERIC registry (`.registry.new`/`register`/`get`/`list`/`unregister`/`dropKind`) + four PER-KIND registries each carrying a versioned CONTRACT — `.model` (pricer), `.signal`, `.calibrator`, `.execution.fillModel` (each with `register`/`get`/`list`/`conforms`). A plug-in registers a HANDLE to the UNCHANGED function + a MANIFEST (declared `in`/`out` schema, dicts of field→typeTag). Conformance (`.registry.conforms` / `.contracts.verify[]`) checks the manifest against the contract (version + every required in/out key present with matching type tag; a superset is fine) and **CAN FAIL** — a deficient/mismatched/wrong-version manifest is rejected (same discipline as the gov gates). Loads right after `config/`, before `models/` (depends on nothing above base q); `core/registry_populate.q` (loaded last) registers the 4 existing capabilities (`.engine.priceOption`/`.strategy.path.commoditySignals`/`.commodity.calibrateCurve`/`.exec.fill`) — METADATA ONLY, no compute change, existing call paths untouched, so the suite stays byte-identical (**383 → 386/0**, +3 registry tests). Strategies keep their own `.strategy.register` (NOT rebased). Contracts documented in `docs/CONTRACTS.md`; demo `apps/examples/registry_inventory.q`. This is the foundation R5 (model cards)/R6 (problem templates)/R7 (agents) plug into. BOUNDED: the mechanism + contracts + registering what exists + a conformance test — NOT a call-site migration.

**Regime layer (`.regime`, v0.64, ARCHITECTURE.md Part II §11.1-11.2 — Research OS R1).** The seam between the compute engine and the research OS: `regime/regime.q` measures a deterministic per-(commodity,date) regime fingerprint from the HDB `futures` data — `curveState` (front-vs-deferred slope), `volState` (realized vol of front-settle log-returns, roll-days zeroed), `liqState` (front volume), `rollPhase` (days-to-front-expiry), `seasonPhase` (reuses `signals/seasonality`) — each with a causal rolling-percentile (`slopePct`/`volPct`/`volumePct`). Measurement only — NO prediction/opinion. API: `.regime.series[commodity;dates]`, `.regime.label[commodity;date]`, `.regime.breakdown[pnlByDate;regimeLabels;axis]` (groups a backtest's daily PnL by regime bucket, **reusing `.strategy.commodityBT.__perf`**, with a `(blended)` row that reproduces the backtest headline). Pure core (`__panelFromLong`/`__labelPanel`) is HDB-free + known-answer tested. Thresholds in `.cfg.regime`. ADDITIVE — depends only on `data/`+`signals/`, never imports `backtest/` (the `__perf` call is lazy); loads after `signals` / before `execution`; never opens the HDB at import; `.regime.open` uses `get` not `\l`. The splayed `regimes` table is built alongside `futures` by `scripts/build_regimes.q` (`.regime.buildTable`). Existing **373 tests byte-identical** (+ 2 synthetic regime tests = **375/0**; canonical numbers unchanged). Demo: `apps/examples/regime_conditional_backtest.q` (the same backtest, now broken down +X in backwardation / −Y in contango). This unblocks R3 (regime-conditional gates), R4 (analogue engine), and the regime library.

**Regime analogue library + risk memory (`.regime.analogue` / `.regime.library`, v0.67, ARCHITECTURE.md §11.3 — Research OS R4).** The "digested history": `regime/analogue.q` adds a curated, versioned library of NAMED historical crude regime EPISODES (the splayed `regimeEpisodes` table + `docs/REGIME_LIBRARY.md`), each with a COMPUTED dominant fingerprint (modal curveState/volState/liqState/rollPhase/seasonPhase + mean slopePct/volPct/volumePct over its date range, from the `regimes` table), its drivers, and its **risk memory** (the known strategy failure modes — "what killed people here"), plus an ANALOGUE engine. `.regime.analogue.distance[a;b]` (pure) = weighted Euclidean over the 3 percentile axes + a per-mismatch categorical penalty over the 5 discrete states (weights in `.cfg.regime.analogue`); `.regime.analogue.nearest[state;n]` / `.regime.analogue.forDate[commodity;date;n]` return the nearest episodes + their risk memory. **Honest data scope:** episodes are matched ONLY on in-HDB windows (crude ≈ 2018.12→2026), each with a real fingerprint; out-of-data lessons (2008, 2014-16) are NARRATIVE ONLY in `REGIME_LIBRARY.md`, never fabricated. **Layer rule:** `analogue.q` reads `regimes`+`regimeEpisodes` and reuses `.regime.label`; it does NOT import `gov/`/`backtest/` (regime stays LOW; gov above may consult the analogue — a legal downward call). Never opens the HDB at import (`.regime.library.open` uses `get`). ADDITIVE: **381 tests byte-identical** (+ 2 synthetic analogue tests = **383/0**). `scripts/build_regime_library.q` builds `regimeEpisodes` (after `build_regimes.q`); demo `apps/examples/regime_analogue_today.q` (latest crude state → nearest episodes + risk memory; today's backwardation/high-vol most resembles the 2022 energy shock, not the 2020 contango crash). Optional later: a gov-side risk-memory check consulting `.regime.analogue` in the thesis/skeptic step.

**Governance layer (`.gov`, v0.65, ARCHITECTURE.md Part II §11.4 — Research OS R3).** The anti-overfitting heart: R1 (`regime/`) SURFACES regime-conditional structure, R3 (`gov/gov.q`) JUDGES it. Four pieces: a hypothesis **registry** (`.gov.register`/`.gov.hypo`), an **append-only trials ledger** (`.gov.logTrial`/`.gov.trials`/`.gov.nTrials` — the honest N, the multiple-testing denominator; a second log ADDS, never overwrites), the **deflated-Sharpe** core (`.gov.psr`/`.gov.deflatedSharpe`, Bailey & López de Prado; `.gov.phiInv` Acklam inverse-normal added, `Phi` **reuses `.validation.__normalCdf`** — the consistency trap is handled: `__perf`'s Sharpe is ANNUALISED, the DSR needs PER-PERIOD = annualised/√annDays), and an ordered **gate cascade** (`.gov.evaluate`, stop-at-first-failure): Gate 0 thesis (+ POST-HOC flag when the tested bucket ∉ claimedRegimes), Gate 1 cost (net Sharpe ≥ hurdle under a REALISTIC `.exec` cost from `.cfg.gov`), Gate 2 deflated Sharpe (DSR ≥ 0.95, N from the ledger — where small samples die), Gate 3 walk-forward (**reuses `.strategy.commodityBT.__splits`** causally as `.alloc.compare` does). Verdicts: `pass`/`reject`/`research`/`regimeConditional`. The non-optional wrapper `.gov.run[hypoId;pnlByDate;regimeLabels;axis]` ALWAYS logs a trial per bucket THEN evaluates — logging can't be silently skipped and **the backtest engine is never edited** (byte-identity). `gov/` is a HIGH layer (MAY import `backtest/`/`regime/` — downward is legal), loads after `backtest/`+`portfolio/` / before `apps/`, never opens the HDB at import (`.gov.open` uses `get` not `\l`; `.gov.flush` persists). ADDITIVE: existing **375 tests byte-identical** (+ 3 synthetic gov tests = **378/0**; canonical numbers unchanged). Thresholds in `.cfg.gov`. Demo: `apps/examples/gov_gate_momentum.q` (real CRUDE momentum, NET of realistic cost — see R3b for the full-cascade verdicts).

**Governance R3b (`.gov` v0.66, ARCHITECTURE.md §11.4 — completes the cascade + fail-safe verdicts).** Adds the sealed-holdout zone, the one-shot holdout gate (Gate 4), and a fail-safe verdict — PURELY ADDITIVE (gov-side date-range slicing only; no engine/backtest/regime edits). (1) **Three data zones** (`.gov.zone.boundaries` pure / `.gov.zone.range[commodity;zone]`): split the sorted dates into train/validate/holdout (`.cfg.gov.zones` = 0.6/0.2/0.2; holdout = most-recent, OOS in time); gates 0-3 + all exploration see `trainValidate` ONLY. (2) **The seal**: `.gov.holdout.read[commodity]` is the SOLE reader of the holdout range. `.gov.holdoutGate[hypoId;runner]` (Gate 4) is ONE-SHOT — scores the strategy on the sealed holdout via `runner[from;to]`, records the look immutably (`holdoutUsedAt` + result on the hypothesis, a ledger row dataZone=`holdout); a second call returns the recorded verdict WITHOUT recomputing/re-reading. (3) **Fail-safe verdict**: a gate FAILURE → `reject` (Gate 0/1) or `research` (Gate 2/3/4) — NEVER tradeable; only an ALL-PASS disposition is tradeable (`pass`, or `regimeConditional` = cleared every gate but edge is regime-gated). A **`tradeable` boolean** (true IFF every gate incl. holdout passed) is added to the verdict; downstream reads it, never the label. FIXED the R3 hole: a deflation-FAILED slice is now `research`, not the tradeable `regimeConditional`. `.gov.runFull[hypoId;runner;axis]` is the complete 5-gate orchestrator: restrict to trainValidate → gates 0-3 per bucket → holdout ONLY if a bucket cleared 0-3 (the seal protects holdout from unearned looks). ADDITIVE: **378 tests byte-identical** (+ 3 synthetic gov tests: zones, holdout one-shot, verdict fail-safe = **381/0**; `test_gov_gates` updated for the corrected mapping). Demo (real CRUDE): backwardation fails cost → `reject`, contango/flat fail deflation → `research`, all `tradeable=0b`, holdout NEVER read (no bucket earned it — the seal held). **Deferred:** the "what's priced in" gate (needs surface/positioning data the HDB lacks — no faked stub).

**Model cards (`.cards`, v0.69, ARCHITECTURE.md knowledge plug-ins — Research OS R5).** The knowledge plug-in that ties the spine together: `cards/cards.q` is a structured, queryable MODEL CARD per capability synthesising what it IS (the R2 contract), what it assumes, which of the three edge sources it claims (`riskPremium`/`structural`/`informational`; `na` for pricers/fill models), where it's valid (R1/R4), and — the honest part — its VALIDATION STATUS **derived from the gov ledger**, never asserted. `.cards.get`/`.list`/`.forCapability` (params are `capName`/`capKind` — NOT the column names, or qSQL `where capabilityName=capabilityName` self-compares); `.cards.validationStatus[capName]` reads the card's `govHypoId` → recomputes the headline deflated Sharpe over the family (reusing `.gov.deflatedSharpe`) + the holdout outcome → `ungated`/`inResearch`/`validated`/`rejected` + a `tradeable` bool (true IFF the sealed holdout passed — can never claim "validated" when the ledger says otherwise); `.cards.riskMemory[capName]` links the R4 regime library. `.cards.audit[]` (pure core `__auditAgainst`) CAN FAIL: it flags undocumented registered capabilities, cards missing a required section, signal/strategy cards with no named edge, and orphan cards — same honesty discipline as the gov gates / R2 conformance. HIGH layer (reads registries + gov + regime library; imported by nothing below; loads after gov, before apps); never opens the HDB at import (`.cards.open` uses `get`). Curated content in `.cfg.cards` (7 cards: the 4 R2 caps + gammaScalp/shortVariance/timeSeriesMomentum), narrative in `docs/MODEL_CARDS.md`, persisted by `scripts/build_model_cards.q`. ADDITIVE: **386 tests byte-identical** (+ 3 cards tests = **389/0**). Demo `apps/examples/model_card_show.q` (one screen: contract + assumptions + edge + live gov-derived validation + risk memory + audit). Foundation R6 (problem-templates) / R7 (agents) read from.

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