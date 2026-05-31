# gov/ — Research governance (ARCHITECTURE.md Part II §11.4, R3 + R3b)

## Purpose
The anti-overfitting heart. R1 (`regime/`) **surfaces** regime-conditional structure honestly; `gov/` **judges** it and refuses to be fooled by a tempting small-sample slice. It provides a hypothesis **registry** (pre-registration), an **append-only trials ledger** (the honest `N` — the multiple-testing denominator), the **deflated-Sharpe** statistic (Bailey & López de Prado), and an ordered, **fail-safe** **gate cascade** (thesis → cost → deflated Sharpe → walk-forward → **sealed holdout**). The cascade is generic: it gates any (hypothesis + net returns + trial history); the crude momentum-by-regime finding is just its first user. **R3b** completed the cascade with three data zones (train / validate / **sealed holdout**), the **one-shot holdout gate** (Gate 4 — each hypothesis sees the holdout exactly once, ever), and a **fail-safe verdict** (`tradeable` is true IFF every gate passed; a gate failure can never be tradeable).

## Dependencies
A **HIGH** layer — dependencies flow downward, so `gov/` **may** import `backtest/` (it sits above it; the bug would be `backtest/` importing `gov/`). It reuses `.strategy.commodityBT.__perf` (Sharpe + the per-period↔annualised conversion) and `.strategy.commodityBT.__splits` (causal walk-forward index splits — the same pattern `.alloc.compare` uses), `.validation.__normalCdf` for `Phi` (the BS/pricing-core normal CDF, **not** reimplemented), and judges the regime-conditional results from `regime/`. Loads after `backtest/` + `portfolio/`, before `apps/`. **Never opens an HDB at import** (`.gov.open` uses `get`, not `\l`); a fresh process gets empty in-memory tables via `.gov.init`.

## Modules
- `gov.q` — `.gov.*`: the PSR/DSR core, the registry + ledger, the data zones, the gate cascade (incl. the one-shot holdout gate), and the `.gov.run` / `.gov.runFull` wrappers.

## Key API
- `.gov.psr[srHat;srStar;n;skew;kurt]` / `.gov.deflatedSharpe[srHat;n;skew;kurt;N;V]` — pure Probabilistic / Deflated Sharpe Ratio. `.gov.phiInv[p]` — inverse normal CDF (Acklam), added here because no inverse-normal existed in the codebase.
- `.gov.register[hypo]` → hypoId (pre-register; idempotent by hypoId). `.gov.hypo[hypoId]` → the registered dict.
- `.gov.logTrial[trial]` → trialId (**APPEND-ONLY** — a second log adds a row, never overwrites). `.gov.trials[hypoId]` / `.gov.nTrials[hypoId]` (= `N` for deflation).
- `.gov.zone.boundaries[dates]` (pure) / `.gov.zone.range[commodity;zone]` → the `[from;to]` for `train`/`validate`/`holdout`/`trainValidate`. Gates 0-3 and all exploration use `trainValidate` ONLY; holdout is the most-recent slice (OOS in time).
- `.gov.evaluate[ev]` → a verdict record (gates 0-3; stops at the first failure; carries `tradeable`).
- `.gov.holdout.read[commodity]` — the **sole** access path to the holdout range. `.gov.holdoutGate[hypoId;runner]` — Gate 4, **one-shot**: scores the strategy on the sealed holdout, records the look immutably; a second call returns the recorded verdict without recomputing.
- `.gov.run[hypoId;pnlByDate;regimeLabels;axis]` → one verdict row per regime bucket (gates 0-3 only, no holdout). **The non-optional path**: it ALWAYS logs a trial per bucket THEN evaluates.
- `.gov.runFull[hypoId;runner;axis]` → **the complete fail-safe cascade**. `runner[from;to]` is the strategy restricted to a date range (the seam that lets gov control which dates the strategy ever sees). Restricts to `trainValidate`, runs gates 0-3 per bucket, and ONLY if a bucket clears 0-3 does the hypothesis earn the one-shot holdout look. Final per-bucket verdict carries `tradeable` (true IFF gates 0-3 cleared AND the holdout passed).
- `.gov.open[hdbPath]` / `.gov.flush[hdbPath]` — load / persist the registry + ledger splays (via `get`/`set` + `.Q.en`).

## The gates (ordered; stop at first failure)
0. **thesis** — registered rationale + a valid `edgeSource` + `claimedRegimes`. Also flags **post-hoc** when the tested bucket is not in `claimedRegimes` (a data-snooping warning that hardens the deflation reading).
1. **cost** — net-of-execution annualised Sharpe ≥ hurdle, after a **realistic** `.exec` cost config (`.cfg.gov.exec`, not the frictionless default).
2. **deflated Sharpe** — `DSR ≥ 0.95`, with `N` = the ledger trial count and `V` = the variance of the family's per-period Sharpes. The bar **rises with N** and **falls with short n** — this is where a tempting small-sample slice dies.
3. **walk-forward** — OOS Sharpe sign-stability across causal folds (reuses `__splits`).
4. **sealed holdout** (R3b) — a ONE-SHOT net-Sharpe test on the most-recent, never-seen-during-development holdout zone. Each hypothesis sees the holdout exactly once; the look is recorded immutably and a second call returns the recorded verdict. Only reached if gates 0-3 cleared.

Verdicts (**fail-safe**): a gate FAILURE → `reject` (Gate 0/1: no mechanism / no edge after cost) or `research` (Gate 2/3/4: plausible thesis, failed evidence) — **never tradeable**. An ALL-PASS disposition → `pass`, or `regimeConditional` (edge concentrated in a non-claimed regime) — the only tradeable verdicts. The verdict record carries a **`tradeable` boolean** (true IFF every gate passed); downstream must read `tradeable`, never infer from the label string.

## Notes
- Thresholds + the realistic cost config live in `.cfg.gov` (no magic numbers). The pure statistical functions are known-answer tested; the gate cascade and ledger are tested on synthetic in-memory tables (`tests/gov/`, no HDB).
- **Consistency trap (handled):** `__perf` reports an ANNUALISED Sharpe; the DSR needs the PER-PERIOD Sharpe with `n` = number of returns, so the layer converts `annualised / sqrt(annDays)`.
- The persisted `hypotheses` / `trials` splays live alongside `futures` / `regimes` under `.cfg.paths.hdb` (gitignored). The append-only contract holds across runs via open (load prior) → logTrial (append) → flush (write back). The holdout one-shot is recorded on the hypothesis (`holdoutUsedAt` + result columns).
- **The seal** is enforced at the API level: `.gov.holdout.read` is the single function that surfaces the holdout range, called only by the holdout gate; everything else restricts to `trainValidate`. The zone discipline is gov-side date-range slicing — the existing data/backtest queries are unchanged (byte-identical).
- ADDITIVE: the engine/backtest/regime code is untouched, so the full suite stays byte-identical. Demo: `apps/examples/gov_gate_momentum.q` (run through `.gov.runFull`).
- **Deferred (later, see ARCHITECTURE.md §13):** the "what's PRICED IN" gate — a faithful version needs options-surface / positioning data the HDB does not have; not faked with a stub.
