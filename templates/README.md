# templates/ — Problem templates: the research-shape plug-in (ARCHITECTURE.md Part II three plug-in kinds, R6)

## Purpose
A **problem template** is a research SHAPE that composes capabilities, declares its own inputs/outputs/validation, and produces a per-day PnL series that flows into the universal gov gates. It is the **third plug-in kind** — CAPABILITY (R2 contracts) · KNOWLEDGE (R5 cards) · **PROBLEM-TEMPLATE (here)** — registered through R2's `.registry.*` as a `template` kind, so a template is a first-class, conformance-checked, card-able plug-in. R6 generalises the directional strategy shape into this abstraction and proves it generalises by shipping one genuinely different shape (relative-value).

## Dependencies
Composes `signals/`/`backtest/`/`execution/` (loads after them). **Independent of `gov/`** — no mutual import: a template *produces* a PnL series + a validation manifest; the gov gates *consume* a PnL series (as in R3). Loads after `backtest/`+`portfolio/`, before `gov/`; the `apps/` demo wires template output → `.gov.runFull`. Never opens the HDB at import.

## Modules
- `template.q` — `.template.*`: the abstraction (the `template` kind + contract on `.registry.*`, `register`/`get`/`list`/`conforms`/`run`) and **template #1 `directionalSignal`** (a faithful wrapper).
- `relative_value.q` — `.template.rv.*`: **template #2 `relativeValue`** — a new shape with its own spread-stationarity gate.

## The template contract + run
A template run returns `pnl`validation`meta: `pnl` is a `(date;pnl)` table (the per-day NET PnL the gov gates judge), `validation` carries the template-specific check results, `meta` the shape + provenance. The `template` contract enforces this uniform **output**; inputs are template-specific (declared in each manifest for documentation, not contract-required). A deficient template manifest (missing `pnl`/`validation`/`meta` in `out`) is rejected by `.contracts.verify` — the contract bites for templates too.

## Template #1 — directionalSignal (the faithful wrapper)
Expresses the existing directional shape as a template by **delegating** to the unchanged signal + backtest functions (`.strategy.path.commoditySignals` + `.strategy.runAndSummarize`); it does **not** reimplement the logic, so its per-day PnL is **byte-identical by construction** to the existing backtest path. Proven in `tests/templates/test_directional_template_faithful.q` — the test that keeps the abstraction honest (it can't silently change the directional numbers). The existing direct call paths stay untouched; the template is an *additional* way to express directional.

## Template #2 — relativeValue (a genuinely different shape)
A two-leg spread mean-reversion shape, generic over a spread definition (the demo instantiates a crude **calendar spread** — front vs a deferred tenor, straight from the curve the regime layer measures). It composes a new mean-reversion signal (fade extremes — the opposite of momentum) and the existing fill model (`.exec.fill`) into a per-day NET PnL series. **Edge source: structural** (spread reversion).
- **Template-specific validation that CAN FAIL:** a spread-stationarity gate (`.template.rv.stationarity` — an AR(1) / OU mean-reversion test: `b<maxAr1` ⇒ mean-reverting, half-life `-log2/log b`). A non-stationary (random-walk, `b≈1`) spread **fails this gate and the RV thesis is rejected before the universal gates even run** — a gate the directional template does not have. This is new code producing new (non-canonical) numbers; it touches no pinned value.

## Notes
- Config in `.cfg.templates.rv` (lookback / entryZ / notional / txnCostRate / deferredIdx / maxAr1). **Not tuned to pass** — whatever verdict the gates return is the system working.
- Demo `apps/examples/relative_value_template.q` (real CRUDE): the crude calendar spread passes its stationarity gate, then flows through `.gov.runFull` (the same cost/deflation/walk-forward/holdout cascade the directional shapes face) for an honest verdict.
- The precondition for vol / optimal-execution / deep-hedging templates later, and the structure R7's agents instantiate when they propose research. (Deep-hedging/control is a named deferral, not built here.)
