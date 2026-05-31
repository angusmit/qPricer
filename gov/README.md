# gov/ — Research governance (ARCHITECTURE.md Part II §11.4, R3)

## Purpose
The anti-overfitting heart. R1 (`regime/`) **surfaces** regime-conditional structure honestly; `gov/` **judges** it and refuses to be fooled by a tempting small-sample slice. It provides a hypothesis **registry** (pre-registration), an **append-only trials ledger** (the honest `N` — the multiple-testing denominator), the **deflated-Sharpe** statistic (Bailey & López de Prado), and an ordered **gate cascade** (thesis → cost → deflated Sharpe → walk-forward). The cascade is generic: it gates any (hypothesis + net returns + trial history); the crude momentum-by-regime finding is just its first user.

## Dependencies
A **HIGH** layer — dependencies flow downward, so `gov/` **may** import `backtest/` (it sits above it; the bug would be `backtest/` importing `gov/`). It reuses `.strategy.commodityBT.__perf` (Sharpe + the per-period↔annualised conversion) and `.strategy.commodityBT.__splits` (causal walk-forward index splits — the same pattern `.alloc.compare` uses), `.validation.__normalCdf` for `Phi` (the BS/pricing-core normal CDF, **not** reimplemented), and judges the regime-conditional results from `regime/`. Loads after `backtest/` + `portfolio/`, before `apps/`. **Never opens an HDB at import** (`.gov.open` uses `get`, not `\l`); a fresh process gets empty in-memory tables via `.gov.init`.

## Modules
- `gov.q` — `.gov.*`: the PSR/DSR core, the registry + ledger, the gate cascade, and the `.gov.run` logging wrapper.

## Key API
- `.gov.psr[srHat;srStar;n;skew;kurt]` / `.gov.deflatedSharpe[srHat;n;skew;kurt;N;V]` — pure Probabilistic / Deflated Sharpe Ratio. `.gov.phiInv[p]` — inverse normal CDF (Acklam), added here because no inverse-normal existed in the codebase.
- `.gov.register[hypo]` → hypoId (pre-register; idempotent by hypoId). `.gov.hypo[hypoId]` → the registered dict.
- `.gov.logTrial[trial]` → trialId (**APPEND-ONLY** — a second log adds a row, never overwrites). `.gov.trials[hypoId]` / `.gov.nTrials[hypoId]` (= `N` for deflation).
- `.gov.evaluate[ev]` → a verdict record (the ordered cascade; stops at the first failure).
- `.gov.run[hypoId;pnlByDate;regimeLabels;axis]` → one verdict row per regime bucket. **The non-optional path**: it ALWAYS logs a trial per bucket THEN evaluates — so logging cannot be silently skipped and the backtest ENGINE stays untouched.
- `.gov.open[hdbPath]` / `.gov.flush[hdbPath]` — load / persist the registry + ledger splays (via `get`/`set` + `.Q.en`).

## The gates (ordered; stop at first failure)
0. **thesis** — registered rationale + a valid `edgeSource` + `claimedRegimes`. Also flags **post-hoc** when the tested bucket is not in `claimedRegimes` (a data-snooping warning that hardens the deflation reading).
1. **cost** — net-of-execution annualised Sharpe ≥ hurdle, after a **realistic** `.exec` cost config (`.cfg.gov.exec`, not the frictionless default).
2. **deflated Sharpe** — `DSR ≥ 0.95`, with `N` = the ledger trial count and `V` = the variance of the family's per-period Sharpes. The bar **rises with N** and **falls with short n** — this is where a tempting small-sample slice dies.
3. **walk-forward** — OOS Sharpe sign-stability across causal folds (reuses `__splits`).

Verdicts: `pass` (all gates) · `reject` (fails thesis or cost) · `research` (fails deflation while claimed, or fails walk-forward) · `regimeConditional` (fails deflation while post-hoc).

## Notes
- Thresholds + the realistic cost config live in `.cfg.gov` (no magic numbers). The pure statistical functions are known-answer tested; the gate cascade and ledger are tested on synthetic in-memory tables (`tests/gov/`, no HDB).
- **Consistency trap (handled):** `__perf` reports an ANNUALISED Sharpe; the DSR needs the PER-PERIOD Sharpe with `n` = number of returns, so the layer converts `annualised / sqrt(annDays)`.
- The persisted `hypotheses` / `trials` splays live alongside `futures` / `regimes` under `.cfg.paths.hdb` (gitignored). The append-only contract holds across runs via open (load prior) → logTrial (append) → flush (write back).
- ADDITIVE: the engine/backtest/regime code is untouched, so the full suite stays byte-identical. Demo: `apps/examples/gov_gate_momentum.q`.
- **Deferred (R3b / later, see ARCHITECTURE.md §13):** the sealed-holdout zone + one-shot holdout gate (R3b — data-zone segregation deserves its own step; the walk-forward gate provides OOS rigor for now); the "what's PRICED IN" gate (needs options-surface / positioning data the HDB does not have — a faithful version is later, not a faked stub).
