# qFDM / qPricer — Architecture & Engineering Design

**Version:** 0.63 &nbsp;|&nbsp; **Tests:** 373 / 0 green &nbsp;|&nbsp; **Loader:** `\l core/init.q`

**What this document is.** The target-state design *and* the incremental build plan. The codebase
migrates toward it behind the test suite — **every step keeps the full suite green and byte-identical;
nothing is a big-bang rewrite.** The canonical pinned numbers must never move unexpectedly:
`gamma_scalp totalPnl=0.6521844`, `short_variance premium=11.97994`, `crude 63.27 -> 57.68`,
`cross-commodity momentum mean +0.7403373 / 0.75 of cells`, `gas carry -0.0243 -> +0.1773`.

**Scope decision (unchanged, deliberate).** qFDM is a **batch research-and-simulation system**, not a
live trading system. Data arrives from Barchart as historical CSVs; there is no real-time feed,
tickerplant, or low-latency path. End-to-end flow:

```
ingest -> store -> price / calibrate -> (NEW) recognise regime -> signal
       -> backtest (simulated execution) -> (NEW) gate -> PnL / Sharpe / verdict
```

Live trading, if ever added, attaches at the data layer and does not change the layers above. The
regime/state machinery below therefore runs **retrospectively over history** for research; a live
"now-cast" version is a future extension at the data layer only.

**This document has two parts:**

* **Part I — The compute engine.** The foundation. *Built* (migration steps 1–5 complete; 6–7 deferred).
* **Part II — The Research Operating System.** The next frontier. *The new design + its roadmap.* This is
  where "how to find real edge and not overfit" is turned into structure.

---
---

# PART I — THE COMPUTE ENGINE (foundation, built)

## 1. Layered architecture — **DONE (v0.56)**

One monorepo organised into layers with a single iron rule:

> **Dependencies flow downward only.** A layer may call layers below it and must **never** call a layer
> above it. A lower layer importing a higher one is a bug, enforced by review.

This is what makes the system scalable, generically extensible, and developable in parallel: once modules
are decoupled behind layer boundaries, separate workstreams (or separate agents) change different layers
without colliding, with the test suite as the merge gate.

| Layer | Responsibility | Status |
|---|---|---|
| `core/` | math, linear algebra, stats, RNG, the loader (`core/init.q`), the config loader (`core/cfg.q`), logging, IPC helpers | built |
| `config/` | `.cfg` value files (`base.q` + optional `{env}.q` overrides) | built |
| `data/` | Barchart parser (standalone) + the splayed HDB (`.data.hdb.*`) | built |
| `models/` | BS/FDM core + all pricers + pricing domain + `assetclass` routing registry | built |
| `calibration/` | iv, surface, objective, calibrate-curve, Kalman MLE, model quality | built |
| `analytics/` | risk / VaR / scenarios / limits / portfolio / reporting / perf | built |
| `signals/` | seasonality (alpha-signal library) | built |
| `execution/` | daily fill / slippage / cost simulation (`.exec`) | built |
| `backtest/` | strategy engine + commodity suite + walk-forward | built |
| `portfolio/` | cross-strategy allocator (`.alloc`) | built |
| `services/` | optional IPC gateway / HDB service / workers | **reserved (deferred)** |
| `scripts/` | `ingest_hdb.q` (+ reserved CI/pipeline) | partial |
| `apps/` | examples + demos | built |
| `regime/` | **(NEW — Part II)** market-state recognition + regime tagging | **to build (R1)** |
| `gov/` | **(NEW — Part II)** hypothesis registry, trials ledger, deflated Sharpe, gate cascade | **built (v0.65, R3)** |
| `agents/` | **(NEW — Part II)** bounded research-agent prompt definitions | **to build (R7)** |
| `tests/` | flat suite, each test starts `\l core/init.q` | built |

The loader (`core/init.q`) loads layers bottom-up in a **fixed, explicit order — there is no
auto-discovery.** When you add a module you must register it in `core/init.q`. It sets `.qfdm.loaded:1b`
and `.qfdm.version`.

## 2. Configuration — `.cfg`, no hardcoding — **DONE (v0.57–58)**

`.cfg` is populated at startup by `core/cfg.q`: load `config/base.q` (every value = the prior hardcoded
literal, types/key-order preserved), then if `QPRICER_ENV` is set, `config/{env}.q` to override a subset.
Unset `QPRICER_ENV` loads base only — the byte-identical default. Domains: `.cfg.fdm`, `.cfg.iv`,
`.cfg.mc`, `.cfg.calib.*`, `.cfg.analytics.*`, `.cfg.paths`, `.cfg.strategy.*` (every per-strategy config).

Best-practice split by data kind: **numeric/system parameters → q dictionaries in a `.q` file**
(native, typed, computable); **tabular reference data** (calendars, specs, holiday tables) → CSV or HDB
keyed tables; **JSON** only at language boundaries (sharing with Python).

## 3. Data & storage — splayed HDB — **DONE (v0.59)**

* **Splayed, UNPARTITIONED** (revised from the original date-partitioned plan): ~48k rows total, so ~1900
  daily partition dirs would be slow + Windows-hostile. `p#` on the `commodity` column; sym columns
  enumerated via `.Q.en`; written with `set`. Revisit date-partitioning only if the data reaches many
  millions of rows.
* Table `futures` (commodity, contractYM, expiry, firstDate, date, OHLC, settle, volume). Tables
  `curves` / `calibrations` / `backtestRuns` remain **future work** (the `regimes` table in §11.1 and the
  `trials` ledger in §11.4 are the first new ones).
* Ingestion (`scripts/ingest_hdb.q` → `.data.hdb.ingest`) **reuses `.parser.futures.loadAll` verbatim**,
  so the store is byte-identical to the parser by construction. Query layer: `.data.hdb.curveAt` /
  `curveHistory` / `dates`. `.data.hdb.open` uses `get` (not `\l`) so it does not change the process
  working directory. Measured ~4.8× faster data-load than CSV parse.
* The HDB dir (`.cfg.paths.hdb`, default `data/hdb`) is **gitignored**, like the CSVs. CSVs are an
  *ingestion source*, not the store.

## 4. Execution simulation — `.exec` — **DONE (v0.60–63)**

The layer that turns *signal returns* into *tradeable PnL* — the realistic upgrade of the old flat cost
rate, generic across every strategy. Order generation (deltas, respecting lots/limits); fill model
(reference ± slippage scaling with order size vs liquidity, spread, impact; configurable convention);
cost model (commissions, fees, bid-ask, financing/carry); realized-PnL attribution (price PnL, cost drag,
financing → Sharpe/Sortino/maxDD/hit-rate/turnover/**cost-adjusted Sharpe**). Every trade on both desks
routes through `.exec.fill`; **frictionless default is byte-identical**, `deltaPV == stepPnl` preserved.
The shared chokepoints `.strategy.__hedgeStep`/`__hedgeInit` (underlying) and `.strategy.__legCost`
(option legs, premium-scaled slippage) mean wiring once covers every hedged strategy. *The gap between
gross and net-of-cost Sharpe is where most apparent edges die; this layer measures it honestly.*

## 5. Portfolio optimizer — `.alloc` — **DONE (v0.62)**

Allocate across strategies from their net-of-execution return panels, evaluate OUT-OF-SAMPLE via causal
walk-forward, and COMPARE methods. `.alloc.weights[returns;method;cfg]`: `equalWeight` (1/N baseline),
`inverseVol`, `minVariance`, `riskParity` (**default** — uses only the covariance, so it doesn't overfit
noisy expected returns), `maxSharpe`, `meanVariance`. Constraints: long-only + fully-invested (capped
simplex projection), weight cap, turnover penalty, covariance shrinkage. `.alloc.compare` estimates
weights on each split's TRAIN columns only and applies them OOS (**causal** — no future return enters a
weight), ranking methods by OOS Sharpe vs the 1/N baseline — the honest "does optimization beat 1/N
out-of-sample?" deliverable. Namespace deliberately distinct from `.portfolio.*` (option-book) and
`.strategy.portfolio.*` (cross-path ensemble).

## 6. Concurrency / IPC services — **DEFERRED**

`peach` with `-s N` slave threads already handles the embarrassingly-parallel work (backtests across
commodity × strategy × split). The IPC service mesh (HDB service + compute worker pool + gateway, via
`hopen`/`hclose`, sync `h"…"`, async `neg[h](…)`, `.z.pg`/`.z.ps`, deferred-async aggregation) solves a
multi-user / multi-core / always-on need that single-user batch research does not have. **Parked until
production onboarding or hardware beyond the single dev laptop.** Over-engineering is the failure mode —
do not stand up a gateway before it is needed.

## 7. CI + scheduled pipeline — **DEFERRED**

Cloud CI needs a licensed q on the runner plus real compute for the full suite + nightly
ingest→calibrate→backtest→report. **Parked until production / adequate hardware.** Locally, the green
suite run by hand is the merge gate. (The nightly pipeline is re-specified in §12 as part of the Research
OS, since it is where the system runs itself.)

## 8. Conventions & invariants (the rules — apply to EVERY task, Part I and Part II)

These are the standing milestone conventions (source of truth: `CLAUDE.md` + `.claude/skills/q-pricing/SKILL.md`):

* **Behaviour must not change in a restructure.** Canonical numbers stay EXACT unless a change is
  explicitly intended; an unexpected move is a bug to fix.
* **Tests use SYNTHETIC data only.** Real CSVs + the HDB are gitignored and touched only by `apps/`
  examples and `scripts/` — never by tests.
* **Shell / process hygiene.** Run every q script as `q file.q </dev/null` *and* end it with `exit 0;`.
  Quick checks foreground; only the full suite + the cross-commodity example may be backgrounded —
  capture to `scratch/` and WAIT (never start the next command while one runs). No orphan q
  (`Get-Process q`). No `rm`/`rmdir`/cleanup (deny-listed). No `mkdir` (kdb+ creates dirs on write).
* **Testing workflow.** During dev run only the affected group (`q tests/run_group.q <group> -q </dev/null`);
  run the full suite (`q tests/run_all_tests.q </dev/null`, ~7–15 min) only as the pre-commit gate.
  Append new test paths to the right `.test.<group>Files` list in `tests/run_all_tests.q`. The runner
  strips top-level `\l core/init.q` lines and `value`s the rest under one already-loaded init; a test
  fails by `'"reason"`.
* **No auto-discovery.** Register every new module in `core/init.q` in dependency order.
* **q source cautions (silent failures).** No bare `/` separator lines (grep `^\s*/\s*$` before the first
  run); write `exp neg 1f` not `exp -1f`; typed numeric vectors take one trailing suffix (`0.25 0.5 1 2 5f`);
  never shadow built-ins (`value`/`count`/`select`/`type`/`string`, two-letter builtins `ss`/`sv`/`vs`,
  `csv` is a constant); inner lambdas can't see outer locals; use `floor` for nearest-rank percentile
  indices (`` `long$3.8 `` rounds half-up). Namespaces explicit/descriptive; private helpers `.__`;
  camelCase; inputs dicts, outputs dicts/tables; per-trade isolation in batch code.
* **Git / docs.** Never stage/commit without explicit approval; every milestone bumps `.qfdm.version` (in
  `core/init.q`), updates `CHANGELOG.md`, and refreshes **every affected README** (touched layers + root).

---
---

# PART II — THE RESEARCH OPERATING SYSTEM (the next frontier)

The engine prices, calibrates, backtests, and costs. The Research OS is the layer that decides **whether
an apparent edge is real**, and the structure that lets the whole thing grow to cover the entire
commodities market without rework. This is where the hard questions — *how do quants find edge, and not
just guess or overfit?* — become architecture.

## 9. The integration principle: one frozen spine, unbounded plug-ins

**Best practice, stated once and obeyed everywhere: build one small generic spine that never changes, and
grow everything else as plug-ins that attach to it through fixed contracts.** The asymmetry is the whole
point — the spine is small and frozen; the plug-ins are infinite.

```
      You (judgment, go/no-go)
              ^
              |  verdicts flow up
   +----------+------------------------------------------+
   |          GENERIC SPINE  (never changes)             |
   |   fixed research pipeline  (state->regime->model    |
   |     ->backtest->gate->verdict)                      |
   |   governance gates  +  trials ledger                |
   |   regime-tagged HDB (three zones, sealed holdout)   |
   +-----------------------+-----------------------------+
                           ^
            register via   |   contracts / registries
                           |
   +-----------------------+-----------------------------+
   |   PLUG-INS  (grow freely — three kinds)             |
   |   1) capability:  pricer / calibrator / signal /    |
   |        fill model / control solver / ML model       |
   |   2) problem template: directional / relative-value |
   |        / vol / optimal-execution / deep-hedging /    |
   |        factor-decomposition (PCA)                   |
   |   3) knowledge: regime-library entries / rationales |
   |        / risk-memory records / model cards          |
   +-----------------------------------------------------+
```

**Nothing you ever add falls outside the three plug-in kinds.** A *capability* implements a component
contract and self-registers (generalising the existing `.strategy.register`). A *problem template* is a
research **shape** that composes capabilities and declares its own inputs/outputs/validation/gates. A
*knowledge* plug-in is data/docs that grow the system's understanding with no code surgery. When you add
deep hedging next year you add a plug-in; you do not reach into the spine. That is what "scalable and
generic across the whole commodities market" means mechanically — and the discipline (§10) is generic too:
every problem template, however exotic the maths, rides the same gates and the same ledger.

## 10. How to get edge, not overfit (the discipline the spine enforces)

**Money in markets comes from exactly three places, and a real edge is always one of them.**

1. **Risk premium** — you're paid to hold a risk others don't want (commodity roll yield in
   backwardation = payment to warehouse price risk; the variance risk premium = payment for selling crash
   insurance). Persistent, large-capacity, but loses precisely in the bad state.
2. **Structural / behavioural dislocation** — someone trades for non-economic reasons (index
   reconstitution, hedging mandate, margin call, calendar-forced roll); you take the other side. Decays as
   it's arbitraged; capacity-limited.
3. **Informational / processing edge** — you forecast or compute better/faster. Rare; erodes.

**The habit: for every apparent edge, name which of the three it is. If you can't, it's noise.**

**The anti-overfitting chain (sacred — built into the spine so every plug-in inherits it):**

```
hypothesis-first (tiny search space)
  -> pre-register the thesis + claimed edge source  (before touching data)
  -> log EVERY trial automatically                  (the honest N)
  -> deflated Sharpe / multiple-testing correction  (the bar rises with N)
  -> walk-forward + parameter PLATEAU (not peak)     (graceful degradation = real)
  -> regime-conditional evaluation                   (worked everywhere, or one lucky regime?)
  -> "what's priced in" check                        (edge lives only in the residual)
  -> sealed holdout, ONE look                        (final, recorded, immutable)
```

Two structural commitments make this honest rather than aspirational. **The ledger logs every run, not
just the survivors** — the count is the denominator the deflated Sharpe needs; without it you compute your
statistics against the one backtest you remember. **The holdout is sealed** — if an agent or an optimiser
can query it, it will eventually overfit to it, so the separation is technical, not polite.

**Agents are bounded** (§11.7): researcher + validator + skeptic + logger, **never an autonomous trader**.
They propose evidence; the gates evaluate; **you** hold the go/no-go. And capital is staged: paper →
small → scaled. The realistic edge envelope for this batch platform is **regime-conditional risk-premium
harvesting and structural / relative-value strategies, validated to institutional standards** — not fast
mispricing or live now-casting (that needs live data + positioning, which the scope deliberately omits).

## 11. The Research OS layers (the new components)

### 11.1 Regime-tagged data — the seam (build first, R1)

The single technical move that unifies the engine and the Research OS: add a **regime fingerprint** for
every (commodity, date) derived from the existing HDB `futures` data — curve state (backwardation /
contango + slope percentile), realized-vol state, liquidity state, roll proximity, seasonal phase. Once
history is labelled, the backtest can report **by regime** (the difference between "0.8 Sharpe" and
"+1.4 in backwardation, negative in deep contango, worst drawdown in crisis vol" — the second is the
tradeable truth), and the analogue engine can find historical neighbours by distance in label space.

### 11.2 Market State Engine — measurement, not prophecy

`.regime.*` is strictly **deterministic measurement**. It computes the state vector; it does **not**
predict or opine. ("Curve is in the 92nd percentile of backwardation, realized vol is high, volume thin."
Never "oil will rise.") Interpretation lives in the agent layer where it gets challenged. Keeping the
engine cold and reproducible is the guardrail against the agent failure mode of telling a clean story.

### 11.3 Regime / Analogue Library + risk memory — explicit, versioned artifacts

A curated, versioned corpus of labelled historical regimes (crude 2008 / 2014–16 / 2020 negative-WTI /
2022 energy shock; gas winter spikes / shoulder-season collapse / the March-April "widow-maker"; power
heatwave spikes / negative prices). Each entry carries date range, curve/vol/liquidity state, drivers,
and — non-negotiably — its **known strategy failure modes** (this is *risk memory*; it must be written
down, not recalled). The analogue query: given a current/hypothetical state, return the nearest historical
regimes and what happened after — so you can answer the only two questions that matter, *this time is
similar to what, and different because what.*

### 11.4 Governance — registry + trials ledger + gate cascade + drift monitor

**Status: largely DONE (v0.65, R3).** `gov/gov.q` (`.gov.*`) ships the registry, the append-only trials
ledger, the deflated-Sharpe core, and the thesis → cost → deflated-Sharpe → walk-forward gate cascade,
with the non-optional `.gov.run` logging wrapper (the backtest engine is never edited). Deferred: the
sealed-holdout zone + holdout gate (→ R3b), the "priced-in" gate (needs surface/positioning data the HDB
lacks), and the drift monitor. The bullets below are the full target design; the *italic deferred* notes
mark what is not yet built.

* **Hypothesis registry** (`.gov` / HDB table `hypotheses`): every idea pre-registered with its economic
  thesis, claimed edge source, instruments, and a-priori parameter ranges — *before* data.
* **Trials ledger** (HDB table `trials`): the engine writes **every** backtest run automatically (the
  honest N), tagging each with hypothesis id, full parameter vector (each knob marked *free*/P&L-optimised
  vs *fixed*/externally-set), data zone + date range, code/version hash, timestamp, and the full metric
  set — including failed and abandoned runs.
* **Gate cascade** (the promotion funnel): thesis → cost (net-of-`.exec`) → deflated Sharpe (penalised by
  the ledger's N) → walk-forward + stress (params/regimes/commodities) → priced-in → sealed holdout (one
  look). Most candidates die at cost and deflation — that is the system working.
* **Sealed holdout zone** *(deferred → R3b)*: a data-segregation discipline over the HDB — train (open) /
  validate (logged) / holdout (sealed, readable only by the final gate, once per strategy, recorded immutably).
* **Drift monitor** *(deferred)*: promoted strategies watched for edge decay; demoted when they fade.

### 11.5 Model Cards

Every model/signal/strategy carries a documented card: inputs, outputs, parameters, assumptions, and —
declared *before* the backtest — the regimes where it should work and where it should fail. A card with a
blank "fails when…" field is a red flag that you have a curve-fit, not a thesis.

### 11.6 Problem-template abstraction — the genericity

Different commodity edges live in different **problem shapes**, and the template is what makes them all fit
one spine:

```
            Problem template (the contract)
            inputs · method · outputs · gates
                          |
        +-----------+-----------+-------------------+
        v           v           v                   v
   Directional   Factors      Optimal           Deep hedging
   signal->PnL   PCA->curve   execution         NN policy->hedge
                 factors      control->schedule
```

A *signal* is "data → position/PnL". *Optimal execution* is "target + state → schedule" (a control
problem, not a PnL backtest). *Deep hedging* is "book + cost model → hedging policy" (a learned policy).
*PCA* is "curve → factors" (an analysis; its gate is factor stability + economic interpretability, not a
Sharpe). All hang off one contract and therefore ride the same spine, gates, ledger, and agents. **Worked
example — adding deep hedging:** write a template declaring inputs (option book, calibrated dynamics, cost
model), method (train an NN hedging policy under frictions), outputs (hedge schedule + PnL distribution),
gates (beats Black-Scholes delta-hedge net of cost? stable across vol regimes? survives holdout?). The NN
is a *capability* plug-in (a Python sidecar) registering through the contract; the spine does not change.

### 11.7 Agent workforce — bounded, with Claude Code orchestration

Six bounded roles, mapped to the discipline: **Rationale/Proposer** (economic mechanism + edge source) →
**History/Regime** (analogues + risk memory) → **Model** (pick template, fill the card) → **Validation**
(controlled tests, regime-conditional, every run logged) → **Skeptic** (attacks throughout: look-ahead,
crowding, post-hoc thesis, parameter fragility, single-regime concentration) → **Ledger** (records
everything). Start as **prompt files** under `agents/` — you do not need an agent platform to begin.
Scale later with Claude Code orchestration (worktrees / Agent Teams / Dynamic Workflows), using the
existing test suite + the `test-runner` subagent as the merge gate, and `CONTRACTS.md` + `CLAUDE.md` +
the q-pricing skill as the shared constitution. **Agents never trade autonomously** — they propose
evidence; gates evaluate; you decide.

## 12. The research workflow (the operational loop)

**The idea lifecycle** (one hypothesis, birth to verdict):

```
Trigger (human idea | regime flag | nightly schedule)
  -> Frame the hypothesis      (Rationale agent: thesis + edge source; register)
  -> Gather evidence           (Regime agent: analogues + risk memory)
  -> Build and card it         (Model agent: pick template, fill model card)
  -> Controlled validation     (Validation agent: walk-forward, regime, cost)   --> Ledger logs EVERY run
  -> Gate cascade              (thesis->cost->deflation->priced-in->holdout)     <-- Skeptic attacks here
  -> Verdict                   (reject | keep-as-research | regime-conditional | promote)
  -> (loop back: drift / iterate)
```

The flow is one-directional and gated — an idea cannot skip from framed to promoted. The most common
*useful* verdict is not a clean promote but "survives cost and walk-forward, only in high-backwardation →
promote as a **regime-conditional** candidate." Promotion then stages paper → small → scaled.

**The nightly pipeline** (the system running itself — same loop on a timer): ingest → rebuild the HDB
**and recompute regime tags** → recalibrate → rerun the suite regime-conditionally → push the overnight
exploration batch through the gate cascade → flag drift → regenerate the `COMMODITY_DESK` report. What
greets you each morning is a **shortlist that already cleared the gates** plus a drift list — the only way
one person scales to a whole market. (Deferred until §7's hardware/CI note is resolved; for now run by hand.)

## 13. Research OS migration roadmap (incremental, behind the suite, byte-identical)

Same discipline as Part I: each step is additive, keeps the suite green and byte-identical, registers in
`core/init.q`, updates the affected READMEs + `CHANGELOG.md` + this doc, and bumps `.qfdm.version`. **Build
order is chosen for lowest risk and highest unblocking, not ambition.**

* **R1 — Regime layer (the seam). ← START HERE.** New `regime/` layer: Market State Engine v1 (§11.2) +
  the `regimes` table (§11.1) + a regime-conditional performance breakdown that reuses
  `.strategy.commodityBT.__perf`. Purely additive — the backtest engine is untouched (the breakdown is a
  composition over daily PnL + labels), so the suite stays **373/0**, plus new synthetic tests (~375).
  Unblocks R3, R4, R5.
* **R2 — Contracts & registries.** Write `CONTRACTS.md`; generalise the existing `.strategy.register`
  pattern to `.model.register` / `.signal.register` / `.calibrator.register` /
  `.execution.fillModel.register`. Formalises the spine's docking interface. Thin, additive.
* **R3 — Governance. ✅ DONE (v0.65).** `gov/gov.q` (`.gov.*`): the `hypotheses` registry + the
  append-only `trials` ledger (the honest N), the **deflated-Sharpe** core (PSR/DSR, Bailey & López de
  Prado; `.gov.phiInv` Acklam inverse-normal added, `Phi` reuses `.validation.__normalCdf`), and the
  ordered **gate cascade** (thesis → cost → deflated Sharpe → walk-forward; stop-at-first-failure, +
  post-hoc flag). The **non-optional** logging wrapper `.gov.run` logs a trial per regime bucket THEN
  evaluates, so the backtest ENGINE is never touched (byte-identical). Reuses `.strategy.commodityBT.__splits`
  (causal walk-forward, as `.alloc.compare` does) + `__perf`. Demo `apps/examples/gov_gate_momentum.q`
  kills the tempting contango/flat momentum slices on cost + deflation + the post-hoc penalty.
  **Deferred to R3b / later (not built):** the **sealed-holdout zone + one-shot holdout gate** → R3b
  (data-zone segregation over the HDB deserves its own bounded step; the walk-forward gate provides the
  OOS rigor for now); the **"what's PRICED IN" gate** → later (a faithful version needs options-surface /
  positioning data the HDB does not have — deliberately not faked with a stub).
* **R4 — Regime / Analogue Library + risk memory.** `docs/REGIME_LIBRARY.md` + backing tables + the
  analogue-query function (state-space distance, nearest-regime retrieval). Knowledge plug-ins.
* **R5 — Model Cards.** `docs/MODEL_CARDS.md` (one card per model/strategy) + a card registry the gates
  read (Gate 0 = "has a card with a populated failure-mode field").
* **R6 — Problem-template abstraction.** Generalise the strategy contract into a problem template; ship
  one new template beyond directional (e.g. a relative-value or a control / deep-hedging stub) to prove
  the abstraction. Capability + template plug-ins.
* **R7 — Agent workforce.** `agents/` prompt files for the six roles; wire to Claude Code orchestration
  later (worktrees / Agent Teams / Dynamic Workflows) with the suite + `test-runner` as the merge gate.

**Deferred (unchanged):** IPC services and cloud CI / scheduled pipeline — until production onboarding or
adequate hardware (§6, §7).

**High-value first slice for the edge-and-honesty goal: R1 + R3** (regime-conditional evaluation +
governance/ledger). That pair is what converts "backtest more strategies" from a liability into something
safe, and it is the part of this whole design that most separates a modelling demo from quant research.

## 14. Conventions (carried from §8; apply to every Research OS step)

Byte-identical discipline · synthetic tests only · additive (don't touch working engine code paths) ·
register new modules in `core/init.q` (no auto-discovery) · shell/process hygiene (`</dev/null`,
`exit 0;`, foreground quick / background+WAIT full, no `mkdir`, no cleanup, no orphan q) · group-run
during dev, full suite as the gate · q source cautions (no bare `/`, no shadowing, `exp neg`, `floor`
percentile, inner lambdas) · update every affected README + `CHANGELOG.md` + this doc + bump
`.qfdm.version` · never stage/commit without approval. Source of truth: `CLAUDE.md` +
`.claude/skills/q-pricing/SKILL.md`.