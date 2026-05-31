# factor/ — Curve factor decomposition (PCA) (Research OS R8)

## Purpose
`.factor.*` decomposes the futures curve into its principal components — **level / slope / curvature** — plus the residuals. It is the **first new analytical capability built on the completed spine** (R1–R7): not new infrastructure, but the first demonstration that a sophisticated research shape plugs into the finished spine (registered via R2, card-able via R5, gated via R3b, regime-aware via the R4 skeptic).

## Dependencies
Reads the HDB curve accessor **downward** (`.data.hdb.curveHistory` — the same curve the regime layer measures, so the factor view is consistent with the regime view). Loads **before `templates/`** (which compose it). Composes/reads only; changes no compute path; never opens the HDB at import.

## Modules
- `factor.q` — `.factor.*`: the PCA capability + the curve-panel accessors.

## Key API
- `.factor.pca[X;k;cfg]` (pure) → `loadings` (M×k, each column a PC), `scores` (T×k), `residuals` (T×M = the panel minus its k-factor reconstruction), `explainedVar` (k-vector), `colMeans`. Deterministic.
- `.factor.curvePanel[ch;nMat]` (pure) — a constant-maturity-by-RANK level panel (the first nMat tenor-ranked contracts per date). `.factor.changePanel[levels]` — per-maturity log-return changes. `.factor.decompose[commodity;dates;cfg]` — the HDB convenience wrapper (curvePanel → changePanel → pca).
- Registered as the R2 `factor` kind: `.factor.register`/`.get`/`.list`/`.conforms`; the `curvePCA` capability passes `.contracts.verify`.

## The PCA (deterministic power iteration)
q has **no built-in eigendecomposition**, so the top-k PCs come from **power iteration + Hotelling deflation** on the covariance matrix. Determinism is pinned so loadings/scores are byte-identical across runs and known-answer testable:
- a **fixed initial vector** (normalised ones), a **fixed tolerance / max-iter** (`.cfg.factor`), and
- a **sign convention** (force a positive loading on the front maturity — eigenvectors are sign-ambiguous, so pinning the sign is mandatory).

`tests/factor/test_factor_pca.q` builds a synthetic 2-factor panel from known loadings + scores and checks recovery (sign-fixed), explained-variance, ~0 residuals, and byte-identical determinism across two runs.

## The factorRelativeValue template (`templates/factor_relative_value.q`)
Trades the mean-reversion of the curve's **cumulative residual** from its k-factor PCA shape (`.template.factorRv.*`), composing `.factor.*` + the RV mean-reversion engine (`.template.rv.*`) + the fill model. Edge source: **structural**. It carries **two** template-specific gates that can fail before the universal gov gates run: **factor stability** (the top loading is stable across the window's two halves) and **residual stationarity** (the traded deviation mean-reverts).

## Notes
- **Do NOT tune k or the lookback to pass the gates** — k-selection and the lookback are exactly where factor PCA overfits, and the deflated Sharpe + walk-forward + holdout exist to punish that. Config in `.cfg.factor`.
- Carded (R5) with real failure modes (factor instability across regimes, k-selection overfitting, lookback sensitivity, the assumption that residual reversion is tradeable rather than microstructure). Demo: `apps/examples/factor_relative_value.q`.
- The template the stochastic-control / deep-hedging work will follow: a new capability + a new problem-template, the same contracts / cards / gates / bounded workflow.
