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

**Status: DONE (v0.57 + v0.58).** A `.cfg` namespace is populated at startup by `core/cfg.q`: load `config/base.q`, then `config/{env}.q` to override, where the environment is selected by the `QPRICER_ENV` variable (`dev` / `test` / `prod`). Unset `QPRICER_ENV` loads base only — the byte-identical default. All tunable library defaults are externalized: core (`.cfg.fdm`/`.cfg.iv`/`.cfg.mc`), calibration (`.cfg.calib.*`), analytics (`.cfg.analytics.*`), paths (`.cfg.paths`), and every backtest strategy + signal config (`.cfg.strategy.*`, v0.58). `models/`/`signals/` needed no sweep (caller-supplied params).

Format, by data kind (this is the best-practice split):

- **System & numeric parameters** (FDM grid sizes, calibration bounds, txn-cost rate, vol target, roll-days-before-expiry, file paths) -> **q dictionaries in a `.q` file**. Native, typed, computable. Preferred over CSV for config.
- **Tabular reference data** (contract calendars, instrument specs, holiday tables) -> **CSV or HDB tables**, loaded as keyed tables. This is where CSV belongs.
- **JSON** (`.j.k` / `.j.j`) -> only at language boundaries (e.g. sharing config with Python).

**Rule:** no magic numbers in `lib`/module code. Every constant resolves through `.cfg`. The migration replaces them module by module behind the tests.

---

## 3. Data & storage — splayed HDB

**Status: done (v0.59).** Historical data moves off loose `Day_*.csv` files onto an on-disk **columnar kdb+ HDB** — the source of end-to-end query speed.

- **Splayed, UNPARTITIONED** (revised from the original date-partitioned plan): the dataset is ~48k rows total (tiny for kdb+), so ~1900 daily partition directories would be slow and Windows-hostile (tens of thousands of small files + antivirus). A single splayed `futures/` directory is far faster. `p#` (parted) attribute on the `commodity` column; sym columns enumerated via `.Q.en`, written with `set`. If the data ever grows to many millions of rows, revisit date-partitioning then.
- Table `futures` (commodity, contractYM, expiry, firstDate, date, OHLC, settle, volume). (`curves` / `calibrations` / `backtestRuns` tables remain future work.)
- Ingestion (`scripts/ingest_hdb.q` → `.data.hdb.ingest`): **reuses `.parser.futures.loadAll` verbatim** (expiry/firstDate derived exactly as the parser), then enumerates + writes the splay. Everything downstream queries the HDB (`.data.hdb.curveAt`/`curveHistory`/`dates`), which is byte-identical to the parser by construction; the parser is now the CSV-read stage feeding ingestion, not re-run per query. `.data.hdb.open` uses `get` (not `\l`) to avoid changing the process working directory. Measured CRUDE data-load: 1884 ms (CSV parse) → 396 ms (HDB) ≈ 4.8× faster.
- The HDB directory (`.cfg.paths.hdb`, default `data/hdb`) is **gitignored**, like the CSVs. CSVs are an *ingestion source*, not the store. `core/init.q` never opens the HDB at import (library-load independence).

---

## 4. (reserved) Curve & surface construction

Lives in `data/`. Today curves are built from real futures (`.parser.*`) and model-implied via calibration. A classical bootstrapped discount curve is **out of scope** for commodities and not required.

---

## 5. Execution simulation — the research centerpiece

**Status: both backtests wired (v0.60 commodity, v0.61 equity); option-leg costs = step 4c.** `execution/execution.q` (`.exec.*`) is a generic daily fill-and-cost layer: `.exec.fill[order;ctx;cfg]` (pure) returns filled quantity (full unless a participation cap binds), an adverse fill price, and cost components (proportional = `rate*|filled|`; slippage = `|filled|*refPrice*bps/1e4`; optional fixed + size/ADV impact). The default config (`.cfg.exec`) is frictionless and reproduces the legacy flat-cost path byte-identically; realism is opt-in via a per-run `exec` sub-config. The commodity backtest (`.strategy.commodityBT.core*`, order = vol-targeted Δposition, refPrice = front price) and the equity engine (the shared hedge helper `.strategy.__hedgeStep`/`__hedgeInit`, order = dollar notional with refPrice=1, preserving the `deltaPV==stepPnl` identity) both route through it and report GROSS vs NET. Matched to daily settle data — a fill-and-cost model, not a tick/LOB simulator. **Step 4c remainder:** per-strategy equity OPTION-LEG costs are not yet routed through `.exec.fill`.

`execution/` is the layer that turns *signal returns* into *tradeable PnL*. It is the realistic upgrade of today's flat cost rate, and it is generic — one execution layer serves every strategy.

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
2. **Config layer.** Add `.cfg`; replace hardcoded constants module by module. *(DONE — v0.57: core/calibration/analytics + paths; v0.58: all backtest strategy + signal configs.)*
3. **HDB.** Stand up the database + an ingestion script; repoint queries to the HDB; keep CSV ingestion as the source. *(DONE — v0.59: splayed (not partitioned) HDB, `.data.hdb.*` + `scripts/ingest_hdb.q`, examples repointed, ~4.8× faster data-load, byte-identical.)*
4. **Execution layer.** Build `execution/`; route the backtest through it instead of the flat cost rate. *(v0.60: commodity BT wired; v0.61: equity engine wired via the shared hedge helper (frictionless default = byte-identical, identity preserved), gross-vs-net reporting. Step 4c remainder: per-strategy equity option-leg costs.)*
5. **Portfolio optimizer.** Add `portfolio/` (allocation across strategies). *(DONE — v0.62: `.alloc.*` — equalWeight/inverseVol/minVariance/riskParity(default)/maxSharpe/meanVariance + long-only/cap/turnover constraints; causal walk-forward `.alloc.compare` ranks methods OOS vs the 1/N baseline.)*
6. **IPC services** *(optional, last)* — gateway + HDB service + workers, only if always-on / multi-core scale is actually needed.
7. **CI + scheduled pipeline.**

The high-value 80% for the research goal is steps 1-4 plus CI.

---

## 10. Conventions

- **During any restructure, behavior must not change.** The only diffs in a move step are file locations, load paths, and the loader; tests stay byte-identical.
- **q-language traps** are recorded in `.claude/skills/q-pricing/SKILL.md`; **workflow rules** (incl. process/shell hygiene: every script ends `exit 0;`, no orphan background shells) in `CLAUDE.md`.
- **Tests use synthetic data only**; real Barchart data is touched only by `apps/` examples and the ingestion pipeline, never by tests.
- **One-directional dependencies** (§1) are enforced by review: a lower layer importing a higher one is a bug.
