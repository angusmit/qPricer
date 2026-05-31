# Model cards — the knowledge plug-in

**Version:** 0.69 (Research OS R5) · **Mechanism:** `cards/cards.q` (`.cards.*`) · **Structured fields:**
`.cfg.cards` → the `modelCards` table (`scripts/build_model_cards.q`) · **Demo:** `apps/examples/model_card_show.q`.

A model card pulls together, per capability, **what it is** (the R2 contract it satisfies), **what it
assumes**, **which edge it claims**, **where it's valid** (R1/R4), and **what the gates actually granted it**
(the R3/R3b ledger) — with an audit (`.cards.audit[]`) that refuses undocumented or edge-less entries. It is
a record you can query that **cannot overstate itself**: `.cards.validationStatus` is **derived from the gov
ledger**, never asserted — `ungated` when nothing is logged, and it can never read `validated` unless the
sealed-holdout gate actually passed.

## The three edge sources
A signal/strategy card must name one (or it's noise): **riskPremium** (compensation for bearing a risk),
**structural** (a persistent market structure/constraint), **informational** (an information/behavioural
edge). Pricers, calibrators and fill models are infrastructure and carry `na`.

## Validation status (derived, not asserted)
`.cards.validationStatus[capabilityName]` reads the card's `govHypoId`, then from the gov ledger: the trial
count `N`, the headline net Sharpe, the **recomputed** deflated Sharpe over the family (reusing
`.gov.deflatedSharpe`), and the holdout outcome from the hypotheses registry. Status: `ungated` (nothing
logged) · `inResearch` (logged, no holdout look yet) · `validated` (holdout passed) · `rejected` (holdout
failed). `tradeable` is true **iff** the sealed holdout passed.

## Audit coverage (honest by design)
`.cards.audit[]` walks **every** registered capability (the four R2 kinds + the 24 strategies) and flags the
undocumented ones, cards missing a required section, signal/strategy cards with no named edge, and orphan
cards. Today **7 capabilities are carded** (below); the remaining strategies are intentionally **not** carded
yet — the audit reports them as a real coverage gap rather than papering over it (curate, don't auto-spam).

---

## Cards (v1)

### `blackScholesFdm` — model / pricer  ·  edge: `na`
Crank-Nicolson / explicit FDM European-vanilla pricer (`.engine.priceOption`), the engine's reference pricer.
**Assumptions:** lognormal diffusion; constant vol/rate over the grid; European exercise. American/barrier +
local-vol combinations are explicitly rejected. **Failure modes:** mispricing under heavy skew/jumps (use the
jump/stochastic-vol models); grid too coarse near the strike/barrier. **Regimes:** any (regime-agnostic).

### `commoditySignalPath` — signal  ·  edge: `riskPremium`  ·  risk memory: `energyShock2022`
Carry / momentum / Kalman / curve-residual signals from a forward-curve history
(`.strategy.path.commoditySignals`). **Assumptions:** ≥3 curve dates; the curve-history schema
(asofDate/tenor/price/contractYM/expiry); the Kalman state estimated **train-only** (causal). **Edge:** trend
+ carry risk premia in the term structure. **Regimes:** strongest in directional/backwardation trends; weak
in range-bound chop. **Risk memory (see REGIME_LIBRARY.md → energyShock2022):** in a supply-shock
backwardation, late momentum chasers round-tripped the reversal — entry timing, not carry, is the killer.

### `schwartz2Curve` — calibrator  ·  edge: `na`
Two-factor Schwartz-Smith forward-curve calibration (`.commodity.calibrateCurve`). **Assumptions:** ≥3
positive-price curve points (lognormal domain); the schwartz2 model family. **Failure modes:** poor fit on a
strongly seasonal gas curve (overlay seasonality first); non-positive settles break the domain. **Regimes:**
any liquid forward curve (a utility).

### `dailyFillCost` — fillModel / execution  ·  edge: `na`
Generic daily fill-and-cost model (`.exec.fill`): participation cap + slippage (bps of notional) + costs.
**Assumptions:** daily settle prices only (no intraday / no limit-order book). **Failure modes:** understates
cost when a position must trade a large fraction of ADV, or in a liquidity vacuum (cf. the 2020 negative-price
print). **Regimes:** any (a model of execution, not of edge).

### `gammaScalp` — strategy  ·  edge: `riskPremium`
Long-gamma delta-hedged position harvesting realized-over-implied vol. **Assumptions:** continuous delta
hedging at the configured band; realized vol exceeds implied; a liquid underlying. **Failure modes:** bleeds
theta in calm/range-bound markets; hedge slippage erodes the edge. **Regimes:** high-realized-vol.

### `shortVariance` — strategy  ·  edge: `riskPremium`
Sells a delta-hedged ATM straddle to harvest the variance risk premium when IV is rich vs a realized-vol
forecast. **Assumptions:** IV richer than the forecast; straddle held to expiry; modest hedge slippage.
**Failure modes:** blows up in vol spikes (2020/2022-type) — the classic short-vol tail. **Regimes:**
calm / mean-reverting vol.

### `timeSeriesMomentum` — strategy  ·  edge: `riskPremium`  ·  risk memory: `energyShock2022`
Long/short the front by the sign of trailing curve momentum (trend following). **Assumptions:** trend
persistence in the front; vol-targeted sizing; turnover survives realistic execution cost. **Failure modes:**
whipsaws in OPEC-cut-defended ranges; round-trips late entries after a supply-shock spike. **Regimes:**
supply-driven backwardation trends. **Governance:** gated under hypothesis `mom_crude_ts` — per R3b the
regime slices fail deflation/cost, so its derived status is `inResearch` / not tradeable, never `validated`.

### `curvePCA` — factor / capability  ·  edge: `na`  (Research OS R8)
PCA of the curve-change panel into level/slope/curvature + residuals (`.factor.pca`), by deterministic power iteration + Hotelling deflation (q has no built-in eig). **Assumptions:** k factors capture the structure; the loadings are stable over the window; deterministic by a fixed init / tol / sign convention. **Failure modes:** **factor instability across regimes** (the loadings rotate when the curve dislocates — see `crude_2020_covid`); **k-selection overfitting** (too many PCs fit noise); **lookback-window sensitivity** (the covariance, hence the loadings, drift with the window). An analytical capability (regime-agnostic itself); the *use* of its residuals for trading is what carries regime risk.

### `factorRelativeValue` — template  ·  edge: `structural`  ·  risk memory: `covid2020`  (Research OS R8)
Fade the curve's cumulative residual from its k-factor PCA shape — factor-structure reversion (`.template.factorRv.*`), composing `curvePCA` + the RV mean-reversion engine + the fill model. **Assumptions:** the factor structure is stable across the window; the traded residual mean-reverts (rather than being microstructure noise); **k and the lookback are NOT tuned to the gates**. **Failure modes:** the factor structure **breaks down when the curve dislocates** (storage saturation 2020 / supply shock 2022) — exactly when a residual-reversion trade is most exposed; the residual may be microstructure (untradeable after costs), which the cost + deflation gates are there to expose. **Template-specific gates** (run before the universal gates): factor stability + residual stationarity. **Governance:** gated under `rv_factor_crude`; whatever the deflation / walk-forward / holdout cascade returns is the honest verdict — k and the lookback are deliberately left untuned so the gates can do their job.

### `marketState` — state / capability  ·  edge: `na`  (Research OS R9)
The point-in-time Market State for `(asOf, commodity)` (`.state.build`), built through the single as-of accessor `.state.asof` — the **only door to history** for foundation code. Returns only data knowable as of `asOf`, logs provenance, and assembles the curve / calendar spreads / tradable universe / feature-placeholder the evidence layer consumes. **Assumptions:** the HDB `expiry` field is correct; futures take the simple `date<=asOf` path (no revisable data yet); the accessor is the *only* door (the existing direct readers — `regime/`, `backtest/` — are **not yet rebased** onto it). **Failure modes:** **revisable data not yet handled** (only the futures path is exercised; the `validDate` revision path is designed but untested on real revisable series); the **provenance log is advisory until R13** consumes it in the evidence audit; the **universe filter trusts the expiry field**; rebasing the existing readers onto the accessor is a deferral, so look-ahead is prevented *for foundation code that uses the door*, not yet enforced everywhere. A data-access capability (regime-agnostic); the foundation chokepoint that makes point-in-time enforceable by construction.

### `curveEngine` — curve / capability  ·  edge: `na`  (Research OS R10)
The curve engine (`.curve.build`): from the as-of data obtained through R9's single door, the clean curve + calendar spreads + the derived features (roll yield / slope / curvature / contango-backwardation classification) + the parallel/slope/butterfly shock operators + immutable daily snapshots. **Assumptions:** reads through R9's as-of door (so every curve is point-in-time by construction); the slope + classification match `regime/`'s convention and threshold (sourced from `.cfg.regime` so they cannot drift); snapshots assume a deterministic curve. **Failure modes:** **thin-liquidity tenors** distorting slope/curvature (the liquidity filter is a blunt volume floor); **classification threshold sensitivity** near the flat boundary (a curve just inside the band flips on a small move); **`regime/` not yet rebased** onto the engine, so two derivations coexist (mitigated by the agreement test, which proves they agree); snapshot **immutability assumes a deterministic curve** (a differing re-snapshot errors, catching non-determinism). A derived-curve capability (regime-agnostic); the curve source the replay engine, season/carry, PnL attribution, and the first real strategy build on.

### `rollEngine` — roll / capability  ·  edge: `na`  (Research OS R11)
The roll-discipline engine (`.roll.active`): which contract is **active** as of a date, decided **from as-of data ONLY** through R9's door, plus roll events and an analytics-only continuous view. Rules: `days_before_expiry` (default, with a multi-day held→next blend window), `volume_switch`, `fixed_calendar`, and `oi_switch` (designed but flagged unavailable). **Assumptions:** the active-contract choice consumes only `.state.asof`'s `date<=asOf` output, so it is point-in-time by construction; the HDB `expiry` field is correct; the roll rule is a config (`.cfg.rolls`) per commodity with a `default`. **Failure modes:** **`oi_switch` is UNAVAILABLE** — the HDB has volume + expiry but **no open interest**, so requesting it errors rather than fabricating OI; the **continuous series is NON-POINT-IN-TIME** (back-adjustment rewrites past levels at every roll) and is for charting / long-run vol ONLY — **never trade or signal off it** (R13 will flag any strategy that does); the **data edge** (last bars are often expiry-day artifacts with no live deferred contract) is handled by dropping those dates from the event/continuous series. The tradable, point-in-time path is the actual contract `.roll.active` names. A roll-discipline capability (regime-agnostic); the contract-selection source the replay engine (R12) and the first real strategy (R16) build on.

### `curveSeasonality` — season / capability  ·  edge: `structural`  (Research OS R14)
Causal curve/spread seasonality (`.season.features`): the same-calendar-month + same-contract-month z-scores, the seasonal factor, the deseasonalised level, and the seasonally-adjusted slope of the front-deferred spread. The same-calendar-month z IS the signal R16's calendar-spread mean-reversion fades. **Assumptions:** the same-month statistic is **CAUSAL** — same-calendar-month observations up to `asOf`, obtained through R9's door — **never full-sample** (a full-sample seasonal stat is look-ahead); a minimum same-month N (`.cfg.season.minObs`) gates a valid z; reuses R10's front/deferred/slope convention. **Failure modes:** **thin same-month history early in the data** (a calendar month has few historical observations → noisy z; mitigated by min-N → null / low-confidence below the threshold); **seasonal patterns break in structural shifts** (the seasonal factor is stale across a regime change); the **causal trailing stats use less data** than a full-sample estimate (the honest cost of being point-in-time). Distinct from the `signals/` seasonality alpha-library. A feature capability; the signal R16 stands on.

### `carryEconomics` — carry / capability  ·  edge: `riskPremium`  (Research OS R14)
Carry / storage economics (`.carry.features`): implied carry (curve-implied roll yield, positive in backwardation), convenience yield, cash-and-carry fair value, the carry signal (roll yield net of the cost of carry), and an inventory-tightness proxy. **Assumptions:** the implied carry comes from the **curve** (we have it — reuses R10's `rollYield`); the convenience yield + cash-and-carry fair value need a rate `r` + storage cost **supplied from `.cfg.carry`** (config, NOT market-observed). **Failure modes:** the **convenience yield + fair value are ASSUMPTION-DEPENDENT** on the assumed `r`+storage (flagged in `assumptionDependent`) — wrong assumptions bias them directly; the **inventory-tightness is a PROXY** from the degree of backwardation (we have **no** real inventory data — it is never fabricated); **cash-and-carry assumes storability** (it is meaningless for non-storable power). A feature capability; the basis for the carry strategies.

### `pnlAttribution` — attribution / capability  ·  edge: `na`  (Research OS R15)
PnL explain + bucketed curve risk (`.attribution.pnl` / `.attribution.risk`): decompose a replay run's realized PnL into **level / slope / curvature / carry / residual** and compute the bucketed curve delta + roll-down / calendar-spread exposure. The quantitative form of §10's "name which of the three edges" — mostly carry ⇒ risk premium, mostly slope/curvature ⇒ structural normalisation, a large residual ⇒ unexplained. **Assumptions:** the level/slope/curvature axes ARE R10's `.curve.shock.parallel|slope|butterfly` operators (the curve shift at the fixed tenor is projected onto them by least squares — not a second basis); carry is the roll-down along the current curve as the held contract ages; the `residual` is the **plug** (`realized − (level+slope+curvature+carry)`), so the five **reconcile to the realized total by construction** (like R13's `pnlTies`). **Failure modes:** the residual lumps together everything the level/slope/curvature/carry basis does not span (a basis-spread or inter-commodity move lands in residual) — a large residual is **unexplained**, which may be noise OR a real unmodelled factor, so interpret with care; it also absorbs costs/financing; the carry attribution depends on R10's `rollYield` definition; R10's slope/curvature shocks pivot on the front, so a front-only position shows ~0 slope/curvature exposure by construction. A decomposition capability (regime-agnostic); it feeds the model card's economic-rationale with numbers (a consumer concern — `cards/` is not rewired here).

### `seasonalCalSpread` — template  ·  edge: `structural`  ·  risk memory: `energyShock2022`  (Research OS R16 — the capstone)
The first real research output (`.template.scs.run`): fade the crude calendar spread against its **seasonally-adjusted same-calendar-month z-score** (R14) — SHORT when unusually wide for the season (z>+entry), LONG when unusually narrow (z<−entry), EXIT on a z-cross of zero or `maxHold`. Trades the **actual legs** (R11/R10), books NET PnL via the fill model, and emits an R12-shape **run record** so the evidence audit (R13) and PnL attribution (R15) can consume it end-to-end. Edge source `structural` (the spread reverts to its seasonal norm). **Assumptions:** the parameters (legs, z-window, entry/maxHold) are **PRE-REGISTERED** in `.cfg.strategy.seasonalCalSpread` and **NOT optimised to the gates** (tuning a calendar spread is exactly where it overfits — the deflation + walk-forward + sealed holdout exist to punish it); the seasonal z is **causal** (same-calendar-month up to `asOf`, R14). **Failure modes:** the spread can **stop mean-reverting** in a regime shift (the template's own stationarity gate, AR(1)/OU, bites — a random-walk spread is rejected before the universal gates); the **seasonal pattern can break** (a structural dislocation — storage saturation / supply shock — where the spread no longer reverts to its seasonal norm); **near-roll liquidity** on the legs; the edge may be **too small after realistic costs**; the parameters are pre-registered, not swept. Run end-to-end through `.workflow.runReplay` (carded → replay → evidence audit → gate cascade → regime skeptic → attribution → a human-escalation packet) — the verdict is HONEST, and a "research, not tradeable" disposition is the architecture working, not a failure.

---

*To add a card: extend `.cfg.cards` (matching the `modelCards` schema; `capabilityName` must match the
registry name), write its narrative here, then confirm `.cards.audit[]` passes it. Never set a validation
status by hand — it is derived from the gov ledger.*
