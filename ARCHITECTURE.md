# qFDM / qPricer — Architecture & Engineering Design

**Version:** 0.74 &nbsp;|&nbsp; **Tests:** 401 / 0 green &nbsp;|&nbsp; **Loader:** `\l core/init.q`

**What this document is.** The target-state design *and* the incremental build plan. The codebase
migrates toward it behind the test suite — **every step keeps the full suite green and byte-identical;
nothing is a big-bang rewrite.** The canonical pinned numbers must never move unexpectedly:
`gamma_scalp totalPnl=0.6521844`, `short_variance premium=11.97994`, `crude 63.27 -> 57.68`,
`cross-commodity momentum mean +0.7403373 / 0.75 of cells`, `gas carry -0.0243 -> +0.1773`.

**Scope decision (unchanged, deliberate).** qFDM is a **batch research-and-simulation system**, not a
live trading system. Data arrives from Barchart as historical CSVs; there is no real-time feed,
tickerplant, or low-latency path. End-to-end flow:

```
ingest -> store -> price / calibrate -> recognise regime -> signal
       -> (NEW) build as-of market state -> replay backtest (simulated execution)
       -> (NEW) evidence audit -> gate -> PnL / Sharpe / verdict
```

Live trading, if ever added, attaches at the data layer and does not change the layers above. The
regime/state machinery below therefore runs **retrospectively over history** for research; a live
"now-cast" version is a future extension at the data layer only.

**This document has two parts:**

* **Part I — The compute engine.** The foundation. *Built* (migration steps 1–5 complete; 6–7 deferred).
* **Part II — The Research Operating System.** The next frontier. *The new design + its roadmap.* This is
  where "how to find real edge and not overfit" is turned into structure. The **judge** (governance,
  gates, deflated Sharpe, sealed holdout, model cards, bounded agents) is **built (R1–R8)**. The
  **evidence layer** it judges — point-in-time market state, a real curve engine, an event-driven replay
  engine, roll discipline, execution realism, an evidence audit, seasonality and carry, PnL explain — is
  the **next frontier (R9–R16, planned)**, because *a rigorous judge is only as honest as the evidence it
  is handed* (§11.8).

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
| `core/` | math, linear algebra, stats, RNG, the loader (`core/init.q`), the config loader (`core/cfg.q`), the registry spine (`core/registry.q`, `.registry`/`.contracts`), logging, IPC helpers | built |
| `config/` | `.cfg` value files (`base.q` + optional `{env}.q` overrides) | built |
| `data/` | Barchart parser (standalone) + the splayed HDB (`.data.hdb.*`) | built |
| `state/` | **(NEW — Part II)** the as-of accessor (the single door to history) + the Market State object (`.state.*`) — point-in-time discipline made enforceable | **built (v0.74, R9)** |
| `curve/` | **(NEW — Part II)** the curve engine (`.curve.*`) — clean per-date curve, spreads, roll yield, slope/curvature, contango/backwardation, curve shocks; immutable daily snapshots | **planned (R10)** |
| `roll/` | **(NEW — Part II)** config-driven roll rules (`.roll.*`) — active-contract mapping + roll events; trade actual contracts | **planned (R11)** |
| `models/` | BS/FDM core + all pricers + pricing domain + `assetclass` routing registry | built |
| `calibration/` | iv, surface, objective, calibrate-curve, Kalman MLE, model quality | built |
| `analytics/` | risk / VaR / scenarios / limits / portfolio / reporting / perf; **(NEW)** bucketed-curve PnL explain (`.risk.attribution`) | built; attribution **planned (R15)** |
| `signals/` | seasonality (alpha-signal library); **(NEW)** commodity curve seasonality + carry feed the Market State | built; `season/` + `carry/` **planned (R14)** |
| `execution/` | daily fill / slippage / cost simulation (`.exec`); **(NEW)** extended for replay realism (participation cap vs as-of volume, roll-window penalty) | built; extension **planned (R12)** |
| `backtest/` | strategy engine + commodity suite + walk-forward; **(NEW)** event-driven **replay mode** alongside research mode, emitting an auditable run record | built; replay mode **planned (R12)** |
| `evidence/` | **(NEW — Part II)** the deterministic evidence audit (`.evidence.*`) — bites on a replay run record before the gates see it | **planned (R13)** |
| `portfolio/` | cross-strategy allocator (`.alloc`) | built |
| `services/` | optional IPC gateway / HDB service / workers | **reserved (deferred)** |
| `scripts/` | `ingest_hdb.q` (+ reserved CI/pipeline) | partial |
| `apps/` | examples + demos | built |
| `regime/` | market-state recognition + regime tagging + analogue library/risk memory | **built (v0.64 R1 / v0.67 R4)** |
| `factor/` | curve factor decomposition (`.factor`) — PCA level/slope/curvature; the first capability on the spine | **built (v0.73, R8)** |
| `templates/` | problem templates (`.template`) — directional (faithful) + relativeValue + factorRelativeValue | **built (v0.70 R6 / v0.73 R8)** |
| `gov/` | hypothesis registry, trials ledger, deflated Sharpe, gate cascade + sealed holdout | **built (v0.65 R3 / v0.66 R3b)** |
| `cards/` | model cards (`.cards`) — knowledge plug-in: contract + edge + gov-derived validation + audit | **built (v0.69, R5)** |
| `workflow/` | the bounded research loop (`.workflow`) — composes gov/cards/template/regime; escalates to the human | **built (v0.72, R7)** |
| `agents/` | bounded research-agent role prompts + `ORCHESTRATION.md` (docs; not loaded) | **built (v0.72, R7)** |
| `tests/` | flat suite, each test starts `\l core/init.q` | built |

The loader (`core/init.q`) loads layers bottom-up in a **fixed, explicit order — there is no
auto-discovery.** When you add a module you must register it in `core/init.q`. It sets `.qfdm.loaded:1b`
and `.qfdm.version`. The new evidence-layer modules load in dependency order: `state/` and `curve/` and
`roll/` low (above `data/`, below the consumers); `evidence/` high (above `backtest/`, below `gov/`).

## 2. Configuration — `.cfg`, no hardcoding — **DONE (v0.57–58)**

`.cfg` is populated at startup by `core/cfg.q`: load `config/base.q` (every value = the prior hardcoded
literal, types/key-order preserved), then if `QPRICER_ENV` is set, `config/{env}.q` to override a subset.
Unset `QPRICER_ENV` loads base only — the byte-identical default. Domains: `.cfg.fdm`, `.cfg.iv`,
`.cfg.mc`, `.cfg.calib.*`, `.cfg.analytics.*`, `.cfg.paths`, `.cfg.strategy.*`, plus the new evidence-layer
domains `.cfg.rolls` (per-commodity roll rules, R11), `.cfg.costs` (per-commodity execution realism, R12),
and `.cfg.state` (universe/liquidity filters, R9).

Best-practice split by data kind: **numeric/system parameters → q dictionaries in a `.q` file**
(native, typed, computable); **tabular reference data** (calendars, specs, holiday tables) → CSV or HDB
keyed tables; **JSON** only at language boundaries (sharing with Python).

## 3. Data & storage — splayed HDB — **DONE (v0.59)**

* **Splayed, UNPARTITIONED** (revised from the original date-partitioned plan): ~48k rows total, so ~1900
  daily partition dirs would be slow + Windows-hostile. `p#` on the `commodity` column; sym columns
  enumerated via `.Q.en`; written with `set`. Revisit date-partitioning only if the data reaches many
  millions of rows.
* Table `futures` (commodity, contractYM, expiry, firstDate, date, OHLC, settle, volume). Tables
  `curves` / `calibrations` / `backtestRuns` remain **future work** (the `regimes` table in §11.1, the
  `trials` ledger in §11.4, and the new `curveSnapshots` / `rollEvents` / `replayRuns` tables in §11.8 are
  the new ones).
* Ingestion (`scripts/ingest_hdb.q` → `.data.hdb.ingest`) **reuses `.parser.futures.loadAll` verbatim**,
  so the store is byte-identical to the parser by construction. Query layer: `.data.hdb.curveAt` /
  `curveHistory` / `dates`. `.data.hdb.open` uses `get` (not `\l`) so it does not change the process
  working directory. Measured ~4.8× faster data-load than CSV parse.
* The HDB dir (`.cfg.paths.hdb`, default `data/hdb`) is **gitignored**, like the CSVs. CSVs are an
  *ingestion source*, not the store.
* **Point-in-time readiness (R9).** The `futures` settle/volume are already point-in-time by `date` (a
  settle on date D was known at end of D), so the as-of accessor over futures needs only a `date<=asOf`
  filter. Revisable data (fundamentals/inventory) would carry `validDate`/`asof`/`source`/`revision`
  columns and take the latest revision `validDate<=asOf` — the accessor is designed revision-aware so it
  is ready, but those columns are not fabricated for data we do not have.

## 4. Execution simulation — `.exec` — **DONE (v0.60–63)**

The layer that turns *signal returns* into *tradeable PnL* — the realistic upgrade of the old flat cost
rate, generic across every strategy. Order generation (deltas, respecting lots/limits); fill model
(reference ± slippage scaling with order size vs liquidity, spread, impact; configurable convention);
cost model (commissions, fees, bid-ask, financing/carry); realized-PnL attribution (price PnL, cost drag,
financing → Sharpe/Sortino/maxDD/hit-rate/turnover/**cost-adjusted Sharpe**). Every trade on both desks
routes through `.exec.fill`; **frictionless default is byte-identical**, `deltaPV == stepPnl` preserved.
The shared chokepoints `.strategy.__hedgeStep`/`__hedgeInit` (underlying) and `.strategy.__legCost`
(option legs, premium-scaled slippage) mean wiring once covers every hedged strategy. *The gap between
gross and net-of-cost Sharpe is where most apparent edges die; this layer measures it honestly.* **R12
extends** `.exec` for replay realism: participation cap against the *as-of* daily volume, and an extra
roll-window liquidity penalty when the roll map (R11) flags a roll — additive, frictionless default still
byte-identical.

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
* **Testing workflow.** During dev run only the affected group via the read-only `test-runner` subagent
  (`q tests/run_group.q <group> -q </dev/null`); run the **full suite from the MAIN THREAD** as the
  pre-commit gate (`q tests/run_all_tests.q </dev/null`, ~7–15 min, backgrounded to `scratch/`, WAIT) — a
  subagent's backgrounded job dies on turn-end. Append new test paths to the right `.test.<group>Files`
  list in `tests/run_all_tests.q`. The runner strips top-level `\l core/init.q` lines and `value`s the
  rest under one already-loaded init; a test fails by `'"reason"`. Set docs to the REAL resulting count;
  do not anchor on a target.
* **No auto-discovery.** Register every new module in `core/init.q` in dependency order.
* **q source cautions (silent failures).** No bare `/` separator lines (grep `^\s*/\s*$` before the first
  run); write `exp neg 1f` not `exp -1f`; typed numeric vectors take one trailing suffix (`0.25 0.5 1 2 5f`);
  never shadow built-ins (`value`/`count`/`select`/`type`/`string`, two-letter builtins `ss`/`sv`/`vs`/`aj`,
  the keyword **`asof`** — name the as-of parameter `asOf`, not `asof`; `csv` is a constant); inner
  lambdas can't see outer locals; use `floor` for nearest-rank percentile indices (`` `long$3.8 `` rounds
  half-up). Namespaces explicit/descriptive; private helpers `.__`; camelCase; inputs dicts, outputs
  dicts/tables; per-trade isolation in batch code.
* **Git / docs.** Never stage/commit without explicit approval; every milestone bumps `.qfdm.version` (in
  `core/init.q`), updates `CHANGELOG.md`, and refreshes **every affected README + this document + CLAUDE.md**.

---
---

# PART II — THE RESEARCH OPERATING SYSTEM (the next frontier)

The engine prices, calibrates, backtests, and costs. The Research OS is the layer that decides **whether
an apparent edge is real**, and the structure that lets the whole thing grow to cover the entire
commodities market without rework. This is where the hard questions — *how do quants find edge, and not
just guess or overfit?* — become architecture. The **judge** half is built (§9–§13, R1–R8); the
**evidence layer** it judges is the next frontier (§11.8, R9–R16).

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
   |     ->replay->evidence audit->gate->verdict)        |
   |   governance gates  +  trials ledger                |
   |   regime-tagged HDB (three zones, sealed holdout)   |
   +-----------------------+-----------------------------+
                           ^
            register via   |   contracts / registries
                           |
   +-----------------------+-----------------------------+
   |   PLUG-INS  (grow freely — three kinds)             |
   |   1) capability:  pricer / calibrator / signal /    |
   |        fill model / curve engine / factor / ML      |
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
*knowledge* plug-in is data/docs that grow the system's understanding with no code surgery. The
evidence-layer work (R9–R16) adds capabilities (the curve engine, roll, the replay mode, the evidence
audit) and templates (the first real strategy) — it grows the spine's plug-ins; it does **not** reorganise
the spine. That is what "scalable and generic across the whole commodities market" means mechanically.

## 10. How to get edge, not overfit (the discipline the spine enforces)

**Money in markets comes from exactly three places, and a real edge is always one of them.**

1. **Risk premium** — you're paid to hold a risk others don't want (commodity roll yield in
   backwardation = payment to warehouse price risk; the variance risk premium = payment for selling crash
   insurance). Persistent, large-capacity, but loses precisely in the bad state.
2. **Structural / behavioural dislocation** — someone trades for non-economic reasons (index
   reconstitution, hedging mandate, margin call, calendar-forced roll); you take the other side. Decays as
   it's arbitraged; capacity-limited.
3. **Informational / processing edge** — you forecast or compute better/faster. Rare; erodes.

**The habit: for every apparent edge, name which of the three it is. If you can't, it's noise.** PnL
explain (§11.8, R15) turns this from prose into a quantitative attribution — the curve/carry/basis/vol
decomposition is the numerical answer to "which of the three."

**The anti-overfitting chain (sacred — built into the spine so every plug-in inherits it):**

```
hypothesis-first (tiny search space)
  -> pre-register the thesis + claimed edge source  (before touching data)
  -> build an as-of market state                    (NEW R9 — point-in-time, single door)
  -> replay through actual contracts/rolls/fills     (NEW R12 — realistic evidence)
  -> EVIDENCE AUDIT, fail-safe                        (NEW R13 — no look-ahead, PnL ties; bites)
  -> log EVERY trial automatically                  (the honest N)
  -> deflated Sharpe / multiple-testing correction  (the bar rises with N)
  -> walk-forward + parameter PLATEAU (not peak)     (graceful degradation = real)
  -> regime-conditional evaluation                   (worked everywhere, or one lucky regime?)
  -> "what's priced in" check                        (edge lives only in the residual; deferred)
  -> sealed holdout, ONE look                        (final, recorded, immutable)
```

Three structural commitments make this honest rather than aspirational. **The ledger logs every run, not
just the survivors** — the count is the denominator the deflated Sharpe needs. **The holdout is sealed** —
if an agent or an optimiser can query it, it will eventually overfit to it, so the separation is technical,
not polite. And **the evidence is audited before the judge sees it** (R13) — a deflated Sharpe on
look-ahead-contaminated PnL is theater on a fabrication; the audit refuses to pass a contaminated run to
the gates at all. *A rigorous judge is only as honest as the evidence it is handed.*

**Agents are bounded** (§11.7): researcher + validator + skeptic + logger + curator + lead, plus (as the
evidence layer grows) a **Data** agent and a **Backtest-QA** agent — **never an autonomous trader**. They
propose evidence and reason about what the deterministic gates cannot encode; the gates evaluate; **you**
hold the go/no-go. Capital is staged: paper → small → scaled. The realistic edge envelope for this batch
platform is **regime-conditional risk-premium harvesting and structural / relative-value strategies,
validated to institutional standards** — not fast mispricing or live now-casting.

## 11. The Research OS layers (the components)

### 11.1 Regime-tagged data — the seam (built, R1)

A **regime fingerprint** for every (commodity, date) derived from the existing HDB `futures` data — curve
state (backwardation / contango + slope percentile), realized-vol state, liquidity state, roll proximity,
seasonal phase. Once history is labelled, the backtest reports **by regime** ("0.8 Sharpe" vs "+1.4 in
backwardation, negative in deep contango, worst drawdown in crisis vol" — the second is the tradeable
truth), and the analogue engine finds historical neighbours by distance in label space.

### 11.2 Market State Engine — measurement, not prophecy (built, R1)

`.regime.*` is strictly **deterministic measurement**. It computes the state vector; it does **not**
predict. Interpretation lives in the agent layer where it gets challenged. *(Note: the regime layer's curve
derivation is generalised by the new curve engine in §11.8/R10 — the regime fingerprint becomes a consumer
of the curve engine rather than a parallel re-derivation, so the two never disagree.)*

### 11.3 Regime / Analogue Library + risk memory (built, R4)

`regime/analogue.q` (`.regime.analogue.*` / `.regime.library.*`): the splayed `regimeEpisodes` table
(named crude episodes with a COMPUTED dominant fingerprint from the `regimes` table) + `docs/REGIME_LIBRARY.md`
(curated drivers + written risk memory), and the analogue engine (`.regime.analogue.distance`
weighted-Euclidean+categorical-penalty / `.nearest` / `.forDate`). Honest data scope: matching runs ONLY
on in-HDB windows with real fingerprints; out-of-data lessons (2008, 2014-16) are narrative-only. `regime/`
stays LOW. The analogue query answers the only two questions that matter: *this time is similar to what, and
different because what.*

### 11.4 Governance — registry + trials ledger + gate cascade (built, R3 + R3b)

`gov/gov.q` (`.gov.*`): the `hypotheses` registry; the append-only `trials` ledger (the honest N); the
deflated-Sharpe core (PSR/DSR, Bailey & López de Prado; `.gov.phiInv` Acklam, `Phi` reuses
`.validation.__normalCdf`); the ordered **gate cascade** thesis → cost → deflated-Sharpe → walk-forward →
**sealed holdout**; the non-optional `.gov.run` logging wrapper + the `.gov.runFull` orchestrator (the
backtest engine is never edited). R3b added the three data zones (`.gov.zone.*`, gates 0-3 see
train+validate only), the one-shot holdout gate (`.gov.holdoutGate` / `.gov.holdout.read` — one look per
hypothesis, recorded immutably), and fail-safe verdicts (a `tradeable` boolean; a gate failure is never
tradeable). Deferred: the "priced-in" gate (needs surface/positioning data the HDB lacks), the drift
monitor.

**Connective wiring (v0.71, cycle-free).** (A) **carded gating** — `.cards.gatedRun` (in `cards/`, ABOVE
gov) refuses to gate a capability whose card has no populated failure-mode field, then delegates to
`.gov.runFull`; (B) the **regime risk-memory skeptic** — `.gov.run`/`.runFull` attach an informational
`riskMemory` annotation read DOWNWARD from `regime/`, changing no pass/fail; (C) the **complete audit** —
`.cards.audit[]` walks `.registry.kinds[]` dynamically. No edge points upward.

**Where R13 attaches.** The new **evidence audit** (§11.8) becomes a precondition the same way carded
gating is: the workflow runs the replay (R12) → calls `.evidence.audit` on the run record → only on PASS
delegates the PnL to `.gov.runFull`. On FAIL the run is rejected and the gates never run. `gov/` is not
modified; the audit lives in its own `evidence/` layer (above `backtest/`, below `gov/`).

### 11.5 Model Cards (built, R5)

`cards/cards.q` (`.cards.*`): a structured card per capability synthesising the R2 contract + assumptions +
the named edge source + regime applicability + the linked R4 risk memory + a VALIDATION STATUS **derived**
from the gov ledger (never asserted; `tradeable` only if the sealed holdout passed). `.cards.audit[]` CAN
FAIL (undocumented capabilities, missing sections, edgeless signal/strategy cards, orphan cards). A card
with a blank "fails when…" field is a red flag that you have a curve-fit, not a thesis.

### 11.6 Problem-template abstraction — the genericity (built, R6 + R8)

Different commodity edges live in different **problem shapes**; the template makes them all fit one spine.
The `template` kind on `.registry.*` (contract = the uniform `pnl`/`validation`/`meta` output; conformance
bites). Built: **directionalSignal** (faithful wrapper, per-day PnL byte-identical), **relativeValue**
(crude calendar-spread mean-reversion with its own stationarity gate), **factorRelativeValue** (R8 — fade
the curve's residual from its k-factor shape, with factor-stability + residual-stationarity gates). The
shape vol / optimal-execution / deep-hedging will follow: a new capability + a new template, the same
contracts, cards, gates, and bounded workflow.

### 11.7 Agent workforce — bounded (built, R7; two roles added with the evidence layer)

Six bounded role prompt docs (`agents/{researcher,validator,skeptic,logger,curator,lead}.md` +
`ORCHESTRATION.md`, not loaded) + the thin executable loop `workflow/workflow.q` (`.workflow.run`) that
composes `.gov`/`.cards`/`.template`/`.regime` and returns a human-escalation packet. BOUNDED by
construction: no path allocates / trades / marks tradeable without all gates; always defers to the human.

**Two roles added as the evidence layer grows** (still bounded, still file-producing, no engine change):

* **Data agent** — enforces as-of discipline. Audits that the single as-of accessor (§11.8/R9) is the only
  door to history; flags any future-revision, future-volume roll, or settlement leakage. Output: a
  data-quality report.
* **Backtest-QA agent** — the judgment layer on top of the deterministic evidence audit (R13). The audit
  bites on what can be written as a rule (no future date, PnL reconciles, roll respected); the agent
  reasons about what cannot (is this roll rule economically sensible? is the universe filter introducing
  survivorship?). Output: a backtest-QA report.

### 11.8 The evidence layer — making the backtest realistic (R9–R16, planned)

**The thesis.** R1–R8 built a rigorous *judge*. A judge is only as honest as the evidence it is handed.
This section is the *evidence-generator*: point-in-time market state, a real curve engine, an event-driven
replay engine, roll discipline, execution realism, an evidence audit, seasonality and carry, and PnL
explain. The two halves meet at the backtester — the gates evaluate a PnL series, and these components make
that series correspond to something tradeable. Each is an additive plug-in on the existing spine, gated and
carded; none reorganises the repo. The optimal-q *design* is sketched below; the literal q lands in the
milestone prompts.

**(a) Point-in-time / as-of discipline + the Market State object — `.state.*` (R9).** Make the as-of rule a
*single chokepoint*, not a discipline scattered across callers. One accessor — `.state.asof[asOf;commodity]`
(parameter `asOf`, never `asof`, the keyword) — is the **only** door to history: it wraps the existing
`.data.hdb.*` query layer (unchanged), returns only rows with `date<=asOf` (and, for revisable series, the
latest revision `validDate<=asOf` — designed-in, futures use the simple path), and logs a *provenance
record* of what it returned. The **Market State** is a dictionary assembled lazily by the accessor: `date`/
`asOf`, `commodity`, `curve` (the existing curve shape), calendar `spreads`, a `features` sub-dict
(placeholder the curve engine + season/carry fill), the tradable `universe` (contracts not expired as-of),
and refs to the cost model + roll calendar. Signals/templates consume the *state*, never raw prices — which
kills a whole class of look-ahead and roll bugs at the type level. R9 introduces the door without rewiring
the existing direct readers (additive, byte-identical); rebasing `regime/` and `backtest/` onto it is a
later deferral.

**(b) The curve engine — `.curve.*` (R10).** `build[asOf;commodity]` composes the accessor → liquidity
filter → roll mapping → curve construction → spreads → derived features (roll yield, slope, curvature,
contango/backwardation classification, the parallel/slope/butterfly shock operators), returning the curve
component of the Market State. Snapshots are written to a splayed `curveSnapshots` table partitioned by
date and **never rewritten** (a new revision is a new row). This generalises the regime layer's slope/
percentile derivation so the two never disagree about the curve.

**(c) Roll discipline — `.roll.*` (R11).** Roll rules live in `.cfg.rolls` (one entry per commodity:
type/days/window/price-source — days-before-expiry / OI-switch / volume-switch / fixed-calendar /
multi-day). `.roll.active[asOf;commodity]` maps to the active contract(s) **deterministically from as-of
data only** and emits `rollEvents` rows when the active contract changes. The backtester trades the
contract the roll map names; the continuous-adjusted series is a *derived view* for analytics only, never
the thing traded — exactly the common commodity-backtest error this prevents.

**(d) The event-driven replay engine + extended execution — `.backtest.replay.*`, `.exec.*` (R12).** The
keystone. The loop is a deterministic **fold over the trading dates** — `(/)` accumulating a state record
`(book; cash; pnl; log)` across dates, each step a pure function `step[state; asOf] -> state'` that calls
the accessor, the signal/template, the (extended) execution model, and the PnL engine. Purity per step is
what makes the run reproducible and auditable: the entire decision trail is the accumulated `log`, written
as the `replayRuns` record (every data access via the accessor, every fill, every roll event, the position
book, the per-step PnL components). Research mode reuses the *same* signal and cost functions vectorised
over the panel, so a replay result can be diffed against the research result. Replay mode is the promotion
gate (the external standard's "gate 8"). Execution extends `.exec` for participation-vs-as-of-volume and a
roll-window penalty; frictionless default stays byte-identical.

**(e) The evidence audit — `.evidence.*` (R13).** A deterministic function `.evidence.audit[replayRun]` →
a pass/fail report, run by the workflow *between the replay and the gates* (a hard precondition, like
carded gating: FAIL → reject, gates never run). It bites on: no datum used outside the accessor's as-of
provenance (subsumes most look-ahead checks once R9 is the only door); no expired/out-of-universe contract
traded; roll rule respected; costs applied; fills present; **position book ties to fills**; **PnL ties to
positions and market moves**; risk report present; and a cross-check that the run's date range stayed
inside train+validate (the holdout *seal* itself stays in governance, §11.4 — the audit verifies, it does
not re-implement). The accessor *prevents* by construction; the audit *verifies* by reading provenance —
a closed loop. The Backtest-QA agent (§11.7) reasons about what the rule-based audit cannot encode.

**(f) Seasonality + carry — `.season.*`, `.carry.*` (R14).** Small capabilities reading the curve engine's
output and returning features merged into the Market State's `features` dict. Seasonality: same-month and
same-contract-month z-scores, seasonal factors, deseasonalised spread, seasonal curve slope (a same-month
z-score is a group-by over historical same-calendar-month spreads from the accessor — better than a naive
all-months z-score for commodities). Carry: implied carry `ln(F/S)/T`, implied convenience yield,
cash-and-carry fair value, a carry signal (realised roll yield − estimated cost-of-carry), an
inventory-tightness proxy. These are the features the first real strategy (R16) is built on, and they force
every strategy to *attribute* its return to carry / inventory / normalisation / supply shock / seasonal
demand. (Complements the existing `signals/` seasonality alpha-library.)

**(g) PnL explain + bucketed curve risk — `.risk.attribution` (R15, extends `analytics/`).** Given a
position book + the curve, return an attribution dict: PnL split into curve-level / slope / curvature /
carry / residual via the curve-shock operators (b), and risk as bucketed deltas per tenor, plus
calendar-spread / basis / roll-down exposures. This is the quantitative form of "name which of the three
edges" (§10) — it feeds the model card's economic-rationale field with numbers, so "explain why this makes
money" becomes an attribution, not prose.

**(h) The first real research output — `factorRelativeValue`'s sibling, seasonally-adjusted calendar-spread
mean-reversion (R16, a `template`).** Crude M1–M3 (or M2–M6) spread, faded against its *seasonally-adjusted*
same-month z-score (R14), traded through actual contracts with roll discipline (R11), filled with realistic
costs (R12), run through replay → evidence audit (R13) → the full gate cascade → the regime skeptic → the
human. Edge: structural / risk-premium. It forces the whole evidence layer to be correct simultaneously,
which is why the external standard recommends it as the first serious strategy and why it is the first time
the foundation is exercised end-to-end on a plausible edge. *Do not tune the spread legs, the z-score
window, or the thresholds to pass — that is where it would overfit, and the deflation/walk-forward/holdout
cascade exists to punish it.*

## 12. The research workflow (the operational loop)

**The idea lifecycle** (one hypothesis, birth to verdict — now with the evidence layer):

```
Trigger (human idea | regime flag | nightly schedule)
  -> Frame the hypothesis      (Rationale/researcher: thesis + edge source; register)
  -> Card check                (carded? else reject — Gate 0 at the cards->gov boundary)
  -> Gather evidence           (Regime/curator: analogues + risk memory)
  -> Build as-of market state  (NEW R9 — Data agent: point-in-time, single door)
  -> Replay through contracts  (NEW R12 — actual contracts/rolls/fills/costs/book)
  -> EVIDENCE AUDIT            (NEW R13 — Backtest-QA + deterministic audit; FAIL -> reject)
  -> Gate cascade              (thesis->cost->deflation->walk-forward->holdout)   <-- Skeptic attacks
  -> Verdict                   (reject | research | regime-conditional | promote) --> Ledger logs EVERY run
  -> Human go/no-go            (paper -> small -> scaled)
```

The flow is one-directional and gated — an idea cannot skip from framed to promoted, and it cannot reach
the gates on contaminated evidence. The most common *useful* verdict is not a clean promote but "survives
cost and walk-forward, only in high-backwardation → promote as a **regime-conditional** candidate." Nothing
marks a strategy tradeable or allocates capital without the human.

**The nightly pipeline** (the system running itself — same loop on a timer): ingest → rebuild the HDB +
recompute regime tags + rebuild curve snapshots → recalibrate → rerun the suite regime-conditionally → push
the overnight exploration batch through replay → evidence audit → the gate cascade → flag drift → regenerate
the `COMMODITY_DESK` report. What greets you each morning is a **shortlist that already cleared the audit
and the gates** plus a drift list. (Deferred until §7's hardware/CI note is resolved.)

## 13. Research OS migration roadmap (incremental, behind the suite, byte-identical)

Same discipline as Part I: each step is additive, keeps the suite green and byte-identical, registers in
`core/init.q`, updates the affected READMEs + `CHANGELOG.md` + this doc + `CLAUDE.md`, and bumps
`.qfdm.version`. **Build order is chosen for lowest risk and highest unblocking, not ambition.**

**Built (R1–R8):**

* **R1 — Regime layer. ✅ DONE (v0.64).** `regime/` Market State Engine v1 + the `regimes` table + a
  regime-conditional performance breakdown; additive, the backtest engine untouched.
* **R2 — Contracts & registries. ✅ DONE (v0.68).** `core/registry.q`: the generic `.registry.*` + four
  per-kind versioned-contract registries + `.contracts.verify` (bites). Registers the four existing
  capabilities as handles to unchanged functions. Metadata only, byte-identical.
* **R3 — Governance. ✅ DONE (v0.65).** `gov/gov.q`: the `hypotheses` registry + the append-only `trials`
  ledger + deflated Sharpe + the ordered gate cascade (thesis → cost → deflated Sharpe → walk-forward) via
  the non-optional `.gov.run`. The backtest engine is never touched.
* **R3b — Sealed holdout + fail-safe verdicts. ✅ DONE (v0.66).** Three data zones, the one-shot holdout
  gate, `.gov.runFull`, and the `tradeable` boolean (a gate failure is never tradeable). The "priced-in"
  gate is deferred (needs surface/positioning data the HDB lacks — not faked).
* **R4 — Regime / Analogue Library + risk memory. ✅ DONE (v0.67).** `docs/REGIME_LIBRARY.md` + the
  `regimeEpisodes` table + the analogue query. In-HDB windows only; `regime/` stays LOW.
* **R5 — Model Cards. ✅ DONE (v0.69).** `cards/cards.q` + `docs/MODEL_CARDS.md` + the `modelCards` table;
  validation status DERIVED from the gov ledger; `.cards.audit[]` bites.
* **R6 — Problem-template abstraction. ✅ DONE (v0.70).** `templates/`: the `template` kind +
  directionalSignal (faithful) + relativeValue (own stationarity gate).
* **R7 — Agent workforce. ✅ DONE (v0.72).** Six bounded role docs + `workflow/workflow.q`; the
  human-escalation packet; bounded by construction.
* **R8 — Factor decomposition (PCA). ✅ DONE (v0.73).** The first shape ON the spine: `factor/factor.q`
  (`.factor.*`, deterministic power-iteration PCA into level/slope/curvature + residuals, pinned sign
  convention) + the `factorRelativeValue` template (own factor-stability + residual-stationarity gates),
  carded + gated. Demo finds no tradeable edge honestly; k/lookback left untuned.

**Planned (R9–R16) — the evidence layer (§11.8):**

* **R9 — Point-in-time market state + as-of discipline. ✅ DONE (v0.74).** New `state/`: the single as-of
  accessor `.state.asof` (the only door to history, wrapping the unchanged `.data.hdb.*` query layer with a
  `date<=asOf` filter + provenance logging; revision-aware by design via `validDate`) and the Market State
  object `.state.build` (a lazily-assembled dict: curve, spreads, features placeholder, as-of universe,
  cost/roll refs). No-look-ahead invariants as enforceable synthetic checks (`.state.invariant.*`, which R13
  reuses). Registered via R2 (`state` kind), carded via R5. Additive — does NOT rewire the existing direct
  readers (byte-identical). Unblocks R10–R16; R9 + R13 close the loop (the door prevents look-ahead by
  construction, the evidence audit proves it by reading the provenance). (Parameter `asOf`, never the keyword
  `asof`.) Demo `apps/examples/market_state_asof.q`.
* **R10 — Curve engine. ← NEXT.** New `curve/`: `.curve.build` (clean per-date curve, spreads, roll yield,
  slope/curvature, contango/backwardation, curve shocks) + immutable `curveSnapshots`. Generalises the
  regime layer's curve derivation. Carded, conformance-checked.
* **R11 — Roll discipline.** New `roll/`: `.cfg.rolls` + `.roll.active` (as-of-only active-contract
  mapping) + `rollEvents`. Trade actual contracts; the continuous series is a derived view only.
* **R12 — Event-driven replay engine + extended execution.** The keystone. `backtest/` gains a replay mode
  (a deterministic fold over dates emitting the auditable `replayRuns` record — book, fills, roll events,
  per-step PnL); `.exec` extends for as-of participation + roll-window penalty. Replay is the promotion
  gate. Frictionless default byte-identical.
* **R13 — Evidence audit.** New `evidence/`: `.evidence.audit[replayRun]` — the deterministic pre-gate
  (no look-ahead, no expired/out-of-universe trade, roll respected, costs applied, book ties to fills, PnL
  ties to positions+moves, date-range containment cross-check). Wired as a hard precondition (FAIL →
  reject, gates never run), exactly like carded gating; `gov/` not modified. The Backtest-QA agent reasons
  about what the rule-based audit cannot encode.
* **R14 — Seasonality + carry.** New `season/` + `carry/`: same-month z-scores + carry signal + implied
  convenience yield, merged into the Market State features. The features R16 is built on.
* **R15 — PnL explain + bucketed curve risk.** `analytics/` gains `.risk.attribution`: PnL split into
  curve/slope/curvature/carry/residual + bucketed deltas per tenor. The quantitative "name which edge."
* **R16 — First real research output: seasonally-adjusted crude calendar-spread mean-reversion.** A new
  `template`, run replay → evidence audit → gate cascade → regime skeptic → human. The first end-to-end
  exercise of the evidence layer on a plausible edge. Do NOT tune to pass.
* **Then:** re-run the R8 `factorRelativeValue` on the realistic replay foundation — its first verdict that
  is judged on audited evidence rather than vectorised PnL.

**Deferred (unchanged):** the "priced-in" gate; the drift monitor; IPC services and cloud CI / scheduled
pipeline (§6, §7); pricing extensions (commodity swaps, Asian / spread options — §9 of the external
standard); dashboards; fundamental-event and volatility-risk-premium strategies (need data we lack); power
specialization (no power data) — added only when a named edge or a real need justifies them.

**High-value first slice: R9 + R12 + R13** (point-in-time state + replay + evidence audit) — that triple is
what converts every gate verdict from rigorous-in-principle to rigorous-in-fact, and it is the part of this
whole design that most separates a modelling demo from a commodities-desk research system.

## 14. Conventions (carried from §8; apply to every Research OS step)

Byte-identical discipline · synthetic tests only · additive (don't touch working engine code paths) ·
register new modules in `core/init.q` (no auto-discovery) · shell/process hygiene (`</dev/null`,
`exit 0;`, group-run via the read-only `test-runner` during dev / full suite from the MAIN THREAD as the
gate, backgrounded + WAIT, no `mkdir`, no cleanup, no orphan q) · don't anchor on a target test count, set
docs to the REAL count · q source cautions (no bare `/`, no shadowing — including the keyword `asof`, use
`asOf`; `exp neg`, `floor` percentile, inner lambdas) · update every affected README + `CHANGELOG.md` +
this doc + `CLAUDE.md` + bump `.qfdm.version` · never stage/commit without approval. Source of truth:
`CLAUDE.md` + `.claude/skills/q-pricing/SKILL.md`.

**Carried-over non-negotiables, evidence-layer additions:**

1. **Point-in-time, enforced at one chokepoint.** The as-of accessor (`.state.asof`) is the only door to
   history for foundation code; nothing reads raw partitions directly.
2. **The gates are only as honest as the evidence.** No PnL reaches the gate cascade without passing the
   deterministic evidence audit (R13). The audit is fail-safe: a contaminated run is rejected, not judged.
3. **Trade actual contracts.** Roll rules are config-driven and as-of-only; the continuous-adjusted series
   is a derived view for analytics, never the thing traded.
4. **Replay is the promotion gate.** A strategy is not viable until it survives event-driven replay, the
   evidence audit, and the full gate cascade — and a human signs off. Capital stages paper → small → scaled.