# qPricer — Architecture & Engineering Design

**Status:** target-state design. The codebase migrates toward this incrementally, behind the test suite (§9). Nothing here is a big-bang rewrite.

**Scope decision (deliberate):** qPricer is a **batch research-and-simulation system**, not a live trading system. Data arrives from Barchart as historical CSVs; there is no real-time feed, tickerplant, or low-latency path. The end-to-end flow is:

```
ingest -> store -> price / calibrate -> signal -> backtest (with simulated execution) -> PnL / Sharpe
```

Dropping live trading is intentional: it removes the entire real-time stack (feed handler, tickerplant, RDB, latency engineering), none of which serves a backtest-driven research goal. If live trading is ever added, it attaches at the data layer and does not change the layers above.

---

## 1. Layered architecture

The system is **one monorepo** organized into layers with a single iron rule:

> **Dependencies flow downward only.** A layer may call the layers below it and must never call a layer above it.

This rule is what makes the system scalable, generically extensible, and — critically — developable in parallel: once modules are decoupled behind their layer boundaries, separate workstreams (or separate agents) can change different layers at once without colliding, with the test suite as the merge gate.

| Layer | Responsibility | Current modules that move here |
|---|---|---|
| `core/` | math, linear algebra, stats, RNG, the config loader, logging, IPC helpers | cholesky / normal generators / optimizer primitives currently inline |
| `config/` | environment configs + tabular reference data (calendars, specs) | (new) |
| `data/` | Barchart ingestion, curve/surface construction, HDB load/query | `parser.q` (`.parser.crude` / `.parser.futures`) |
| `models/` | pricers: valuation + greeks | `schwartz.q`, `schwartz2.q`, `meanRevertingJump.q`, `commoditySpread.q`, the BS/FDM core |
| `calibration/` | fit model parameters to market; the generic optimizer | `calibrateCurve.q`, `kalmanSchwartzSmith.q`, `calibration.q`, `objective.q` |
| `analytics/` | greeks, scenarios, VaR, P&L attribution, model-risk limits | risk / scenario / limit modules |
| `signals/` | alpha signal library | signal half of `commodityStrategies.q`, `seasonality.q` |
| `execution/` | **(new)** order generation + fill / slippage / cost simulation | (new — see §5) |
| `backtest/` | strategy engine, registry, walk-forward, scoring | `strategy.q`, the strategy definitions |
| `portfolio/` | **(new)** allocation / portfolio optimizer | (new) |
| `services/` | **(optional)** IPC processes: gateway, HDB service, compute workers | (new — see §7) |
| `apps/` | examples, report generators, CLIs | `examples/` |
| `tests/` | mirrors the tree | `tests/` |
| `scripts/` | build, CI, the scheduled pipeline | (new) |

The loader (`core/init.q`, evolving from today's `lib/init.q`) loads layers bottom-up in a fixed, explicit order. Each layer may have a small `_load.q` manifest listing its files.

---

## 2. Configuration — no hardcoding

**Status: in progress (step 2 done; backtest layer deferred to step 2b).** A `.cfg` namespace is populated at startup by `core/cfg.q`: load `config/base.q`, then `config/{env}.q` to override, where the environment is selected by the `QPRICER_ENV` variable (`dev` / `test` / `prod`). Unset `QPRICER_ENV` loads base only — the byte-identical default. As of v0.57 the core (fdm/iv/mc), calibration (curve/Kalman) and analytics (model-report thresholds) defaults plus `.cfg.paths` are externalized; the `backtest/` per-strategy configs remain (step 2b).

Format, by data kind (this is the best-practice split):

- **System & numeric parameters** (FDM grid sizes, calibration bounds, txn-cost rate, vol target, roll-days-before-expiry, file paths) -> **q dictionaries in a `.q` file**. Native, typed, computable. Preferred over CSV for config.
- **Tabular reference data** (contract calendars, instrument specs, holiday tables) -> **CSV or HDB tables**, loaded as keyed tables. This is where CSV belongs.
- **JSON** (`.j.k` / `.j.j`) -> only at language boundaries (e.g. sharing config with Python).

**Rule:** no magic numbers in `lib`/module code. Every constant resolves through `.cfg`. The migration replaces them module by module behind the tests.

---

## 3. Data & storage — partitioned HDB

Historical data moves off loose `Day_*.csv` files onto a **date-partitioned, splayed kdb+ HDB** — the core reason to use kdb+ at all, and the source of end-to-end query speed.

- Partition by date; splay columns; `p#` (parted) attribute on the contract/`sym` column; `g#` where useful.
- Core tables: `futures` (date, commodity, contractYM, OHLC, settle, volume), `curves` (date, commodity, tenor, price), `calibrations` (date, commodity, model, params...), `backtestRuns` (runId, strategy, metrics...).
- Ingestion: parse a Barchart CSV -> write into the HDB via `.Q.dpft`. Everything downstream queries the HDB (columnar, memory-mapped), never re-parses CSVs.
- The HDB directory is **gitignored**, like the CSVs. CSVs are an *ingestion source*, not the store.

---

## 4. (reserved) Curve & surface construction

Lives in `data/`. Today curves are built from real futures (`.parser.*`) and model-implied via calibration. A classical bootstrapped discount curve is **out of scope** for commodities and not required.

---

## 5. Execution simulation — the research centerpiece

`execution/` is the new layer that turns *signal returns* into *tradeable PnL*. It is the realistic upgrade of today's flat cost rate, and it is generic — one execution layer serves every strategy.

1. **Order generation** — target positions -> orders (deltas to trade), respecting lot sizes and position limits.
2. **Fill model** — execution price = reference +/- slippage, where slippage scales with order size relative to liquidity (volume / ADV), spread, and a market-impact coefficient. Configurable fill convention (next-bar open, VWAP, close +/- half-spread).
3. **Cost model** — commissions, fees, bid-ask, and financing / carry on held positions.
4. **Realized PnL & attribution** — decompose realized PnL into price PnL, cost drag, and financing; then compute Sharpe, Sortino, max drawdown, hit rate, turnover, and **cost-adjusted** Sharpe.

The gap between gross signal Sharpe and net-of-cost Sharpe is where most apparent edges die; this layer measures it honestly.

---

## 6. Extensibility — registries & contracts

Generalize the existing `.strategy.register` pattern across the system: `.model.register`, `.signal.register`, `.calibrator.register`, `.execution.fillModel.register`.

- Each component type implements a documented **contract** (a fixed set of functions): a pricer implements `price` / `greeks`; a signal implements `compute`; a fill model implements `fill`; a strategy implements `init` / `step` / `summary` (already true).
- Engines dispatch by registry lookup. A new model/signal/strategy is a new file that **self-registers** — the core never changes.
- Contracts are documented in `CONTRACTS.md`. This is what lets the system grow additively ("develop on its own") rather than invasively.

---

## 7. Concurrency

**Execution parallelism (do this first):** `peach` with `-s N` slave threads for the embarrassingly-parallel work (backtests across commodity x strategy x split). For true multi-core beyond slaves, multiple worker processes coordinated by a gateway via async IPC.

**IPC services (optional, deferred until needed):**

- **HDB service** — keeps the partitioned database memory-mapped and always-on, so it is never reloaded.
- **Compute worker pool** — N stateless q processes; the gateway hands out work (q primitives are single-threaded per process, so multiple processes give real multi-core parallelism beyond `-s`).
- **Gateway** — single entry point; fans work out to workers (async), aggregates, queries the HDB service. Clients talk only to the gateway.
- Mechanics: `hopen` / `hclose`; sync `h"expr"` for queries; async `neg[h](func;args)` for dispatch; `.z.pg` / `.z.ps` handlers; deferred-async with callbacks for non-blocking aggregation in the gateway.

**Honest guidance:** for a single user, a partitioned HDB + `peach` + `-s` slaves in one process already delivers most of the speed. Build the service mesh only when reload time or a single-process compute ceiling becomes a real pain. Do not stand up a gateway before it is needed — over-engineering is the failure mode.

---

## 8. Quality & self-running workflows

- **CI:** run `tests/run_all_tests.q` on every commit / PR (GitHub Actions with a q Docker image, or a local pre-commit hook). Wire in the bare-`/` scanner, a no-orphan-shell check, and formatting checks as gates. The full suite is the guard that makes refactoring safe.
- **Scheduled pipeline (the system maintaining itself):** a cron / q-timer job that nightly ingests new Barchart CSVs -> updates the HDB -> recalibrates -> reruns the strategy suite + walk-forward -> regenerates the `COMMODITY_DESK` report -> flags drift (a strategy's Sharpe moving, a calibration failing).
- **Versioning:** keep `.qfdm.version`; tag releases in git; maintain the changelog in `CLAUDE.md`.

**Supporting tech, kept lean:** q for compute and storage; a thin layer of **Python** for orchestration / scheduling / plotting (already used for bc-utils); **cron / q timer** for automation; q-native config. Resist additional languages or frameworks — speed comes from the HDB and vectorized q, not infrastructure.

---

## 9. Migration roadmap (incremental, behind the test suite)

The ~360-test green suite is what makes a large refactor safe. Each step keeps every test green and byte-identical; do not big-bang.

1. **Layer the folders.** Move files into the layer tree; fix load paths; restructure the loader. Pure move — **no behavior change**. Lowest risk; do first.
2. **Config layer.** Add `.cfg`; replace hardcoded constants module by module. *(v0.57: done for core/calibration/analytics + paths; backtest per-strategy configs deferred to step 2b.)*
3. **HDB.** Stand up the partitioned database + an ingestion script; repoint queries to the HDB; keep CSV ingestion as the source.
4. **Execution layer.** Build `execution/`; route the backtest through it instead of the flat cost rate.
5. **Portfolio optimizer.** Add `portfolio/` (allocation across strategies).
6. **IPC services** *(optional, last)* — gateway + HDB service + workers, only if always-on / multi-core scale is actually needed.
7. **CI + scheduled pipeline.**

The high-value 80% for the research goal is steps 1-4 plus CI.

---

## 10. Conventions

- **During any restructure, behavior must not change.** The only diffs in a move step are file locations, load paths, and the loader; tests stay byte-identical.
- **q-language traps** are recorded in `.claude/skills/q-pricing/SKILL.md`; **workflow rules** (incl. process/shell hygiene: every script ends `exit 0;`, no orphan background shells) in `CLAUDE.md`.
- **Tests use synthetic data only**; real Barchart data is touched only by `apps/` examples and the ingestion pipeline, never by tests.
- **One-directional dependencies** (§1) are enforced by review: a lower layer importing a higher one is a bug.
