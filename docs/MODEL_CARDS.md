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

---

*To add a card: extend `.cfg.cards` (matching the `modelCards` schema; `capabilityName` must match the
registry name), write its narrative here, then confirm `.cards.audit[]` passes it. Never set a validation
status by hand — it is derived from the gov ledger.*
