# qFDM: From Equity FDM to a Commodity Desk

*A narrative walkthrough of the commodity stack, with real numbers harvested from the
example scripts. Every figure below is actual output, not illustrative.*

---

## 1. What qFDM is

qFDM is a kdb+/q pricing and risk framework. Its **equity finite-difference core** —
Black-Scholes, explicit and Crank-Nicolson solvers, American early exercise, knock-out
barriers, local volatility, Greeks, scenario risk — exists as the **correctness-validation
case**: equity vanillas have closed forms and parities, so they prove the numerical
machinery is right. The **target asset class is commodity** (oil today; power/gas next):
forward-curve-based pricing, stochastic convenience yield, mean reversion, and the spikes
that make energy interesting. AAPL option data is a test fixture, not the destination.

## 2. The commodity stack

On top of the equity core sits a commodity term-structure layer with four models —
**Black-76** (futures options), **one-factor Schwartz** (log mean reversion), **two-factor
Schwartz** (a short-term deviation plus a stochastic equilibrium level, i.e. a stochastic
convenience yield), and a **mean-reverting-jump** model for power spikes — and a **generic
strategy engine**. The engine is registry-based: a strategy supplies `init`/`step`/
`summary`/`defaultConfig`, self-registers, and is driven by one untouched fold over the
path. Its defining discipline is the **portfolio-value identity**: every strategy's reported
`stepPnl` is checked by independently rebuilding `PV = cash + legMarks + hedge·spot` from
state components and asserting `ΔPV == stepPnl` to machine epsilon. Accounting is *verified*,
never just *asserted*.

## 3. The spread complex: the trades that make a commodity book

A commodity desk is largely a book of **spreads** — spark (power vs gas), crack (product vs
crude), calendar (near vs deferred). The spread-option pricers and their cross-checks:

| Method | What it is | Validation |
|---|---|---|
| **Margrabe** | Exact exchange option `max(F1−F2, 0)` (zero-strike spread) | closed form |
| **Kirk** | Fast approximation for non-zero strike | reduces to Margrabe at K=0 |
| **MC** | Two correlated lognormal forwards, terminal payoff | unbiased reference |

The three pricers police each other (params F1=75, F2=70, T=0.25, σ₁=0.30, σ₂=0.25, ρ=0.8,
r=0.05):

- **Margrabe call = 5.761143, put = 0.823.** Kirk at strike 0 reproduces the Margrabe call
  with **difference 0.0** (Kirk collapses *exactly* to Margrabe there — not "close", identical).
- **Put–call forward parity** `call − put = e^{−rT}(F1−F2)` holds to a **residual of 7.1e-15**
  (machine epsilon).
- For a strike-3 spread call, **Kirk = 3.727965** and **Monte Carlo = 3.733962** (standard
  error 0.0196): the approximation sits **within ~0.3 standard errors** of the unbiased MC.

The `sparkSpread` / `crackSpread` strategies are long these Kirk-priced spread options on a
**Cholesky-correlated multi-commodity curve**, delta-hedged in both legs under the futures
mark-to-market convention — and they too pass the PV-identity accounting test to ~1e-15.

## 4. Real data: the Barchart WTI parser

The first real-data path ingests Barchart WTI daily settle CSVs (`Day_CRUDE_YYYYMM00.csv`,
one file per delivery month). Parsing decisions worth noting: the date is the **first 10
chars of the ISO `Time`** field (the DST-shifting offset is ignored), the settle is the
`Latest` column, and a contract's **expiry is derived from the data — the max date in its
own file**, not a hardcoded CME rule (verified: the Jan-2020 file ends 2019-12-19, the actual
last trade). The loader ingests **24,364 daily rows across 96 contracts**, spanning
2018-12-27 to 2026-05-29.

`curveAt[asof]` assembles a `(tenor, price)` snapshot from the contracts alive and quoted on
a date. The **WTI forward curve as of 2020-01-06**:

| tenor (yr) | price | contract | expiry |
|---|---|---|---|
| 0.041 | 63.27 | 202002 | 2020.01.21 |
| 0.123 | 63.04 | 202003 | 2020.02.20 |
| 0.203 | 62.69 | 202004 | 2020.03.20 |
| … | … | … | … |
| 0.959 | 57.68 | 202101 | 2020.12.21 |

Front **63.27** declining to **57.68** a year out: **backwardation**. As an external reality
check, WTI did trade around **$63** in early January 2020 — the parser reproduces the actual
market, not a synthetic level.

## 5. Single-curve calibration

**Economically:** the slope of the forward curve encodes the *convenience yield* — the
benefit of holding physical barrels. Backwardation (front above deferred) means the
convenience yield exceeds the cost of carry. **Mechanically:** we fit the two-factor Schwartz
model to that real curve. The trick that makes it robust is that `log F(τ)` is **linear in
the level parameters (X0, Y0, μY) for a fixed κ**, so those are solved by *exact* least
squares (via `inv`), leaving only the single nonlinear parameter κ to a 1-D grid — no fragile
high-dimensional search.

Fitting the Jan-2020 curve (vols fixed at 0.35 / 0.20 / 0.30):

```
shortFactor0 (X0)  = -0.00325      fitRmse = 0.0738   (≈0.12% on ~$60)
longFactor0  (Y0)  =  4.153        riskFreeRate        = 0.02
meanReversionSpeed =  3.0          impliedFrontSlope   = -0.0603
longDrift    (μY)  = -0.142        netConvenienceYield =  0.0803
```

The **economic read is the real check**: recovered **convenience yield 0.080 > rate 0.02**,
consistent with the observed backwardation — and a contango curve flips the sign. RMSE alone
would not tell you the fit is economically sensible; the convenience-yield comparison does.

**The honest caveat:** a single, nearly log-linear backwardated curve **under-identifies κ**.
The fit railed κ to the upper bound (3.0) — the curve's gentle shape is explained almost as
well by the linear drift as by mean reversion, so the snapshot cannot pin the reversion
*speed*. That limitation is exactly what motivates the next section.

## 6. Kalman state-space estimation: identifying κ from the dynamics

A single snapshot sees the curve's *shape*; the **panel of curves over time** sees its
*dynamics* — and the autocorrelation of the short-term factor is what pins κ. We cast the
Schwartz-Smith two-factor model (`ln F = e^{−κ(T−t)}·χ_t + ξ_t + A(T−t)`) as a
**linear-Gaussian state space** and estimate all parameters by **Kalman-filter maximum
likelihood** over the futures panel. The filter is the textbook predict/update recursion,
threaded by a `scan`, with the log-determinant via Cholesky; the MLE is coordinate descent
on the negative log-likelihood, reusing the repo's grid-search paradigm.

**The identifiability proof is the synthetic round-trip:** simulate the model with *known*
parameters, generate a ~104-date panel, estimate, and recover them.

| parameter | true | recovered |
|---|---|---|
| κ (mean reversion) | 1.5 | **1.487** |
| σ_χ (short vol) | 0.25 | 0.258 |
| σ_ξ (long vol) | 0.12 | 0.109 |

κ comes back to within ~1%, **interior to the search box** — the panel identifies what the
snapshot could not. Running the same estimator on the **real WTI 2020-2021 panel** (25 dates,
300 observations):

```
kappa (mean reversion) = 1.9508     rho       = 0.664
sigChi (short vol)     = 0.3974     muXi      = -0.040
sigXi  (long vol)      = 0.2674     measSigma = 0.00658
logLik                 = 937.03
```

κ ≈ **1.95, identified and interior** — not pinned to a bound. The single-snapshot fit of
§5 *could not produce this*; the time series can.

## 7. The 2020-2021 story, told through the model

The convenience-yield time series (calibrated per as-of date across the curve history) reads
the COVID oil shock as a regime narrative. Selected points:

| as-of | front | convenience yield | regime |
|---|---|---|---|
| 2020-01-02 | 61.18 | +0.092 | backwardation |
| 2020-02-03 | 50.11 | −0.019 | contango |
| **2020-04-02** | **25.32** | **−0.674** | contango |
| **2020-05-04** | **20.39** | **−0.959** | contango |
| 2020-12-31 | 48.52 | +0.001 | contango |
| 2021-02-02 | 54.76 | +0.087 | backwardation |
| 2021-07-02 | 75.16 | +0.184 | backwardation |

Backwardation in January (convenience yield > carry) flips to contango in February as demand
collapses, deepens to **super-contango** at the April–May storage crisis (a deeply negative
convenience yield of **−0.67 then −0.96**, with the front at **$25 then $20** — the market
paying enormously to defer delivery because storage was full), then recovers back to
backwardation by February 2021. The two **regime transitions** are dated automatically:

```
2020.02.03   backwardation -> contango     (COVID demand collapse)
2021.02.02   contango      -> backwardation (recovery)
```

The Kalman **filtered factors** tell the same story structurally: the short-term deviation χ
plunges to **−0.70 in May-2020** (the transient shock) while the equilibrium level ξ falls
from **4.17 (≈$64) to 3.69 (≈$40)** and recovers to **4.30 (≈$74)** — a permanent repricing
distinct from the transient dislocation. The model *decomposes* the shock into "temporary
panic" (χ) and "the world has changed" (ξ).

The deeply negative **April-20-21 settlement prices are excluded** throughout: the Schwartz
models are log-price and cannot represent a non-positive price, so those dates are guarded out
of the lognormal domain rather than fed to the fit.

## 8. Does it trade? — out-of-sample strategy verdict

The point of the models is signals. Eight strategies trade the WTI curve on a
**causal signal-augmented path**: a roll-adjusted continuous front-return series
plus per-date convenience-yield, Kalman χ/ξ, momentum and curve-residual signals
— each computed using only data up to that date. The Kalman parameters are
estimated **once on a training window** and filtered forward; positions are
**vol-targeted** to a common 15% (so the ranking compares edge, not leverage),
**lagged one step**, charged **transaction costs**, and booked under the
futures mark-to-market identity (`PV = cash`, `ΔPV = stepPnl`).

**A single split flatters to deceive.** On the obvious split — train 2020, test
2021 — the ranking is:

| strategy | 2021 OOS Sharpe |
|---|---|
| curveRelativeValue | **2.19** |
| timeSeriesMomentum | 0.48 |
| twoTimescale | 0.42 |
| carryMomentumCombo | 0.12 |
| convenienceYieldCarry | −0.31 |
| chiReversion | 0.0 (didn't trigger) |
| storageCashCarry | −2.63 |

Taken alone this says "relative value wins." **It doesn't survive walk-forward.**
Re-running over **5 rolling splits** (180-trading-day train, 63-day test, sliding
across 2020-01 → 2021-12), the distribution of out-of-sample Sharpe is:

| strategy | mean OOS Sharpe | dispersion (std) | splits positive |
|---|---|---|---|
| **twoTimescale** | **+1.77** | 1.28 | **4 / 5** |
| chiReversion | +0.94 | 0.63 | 4 / 5 |
| carryMomentumCombo | +0.44 | 1.75 | 3 / 5 |
| timeSeriesMomentum | +0.43 | 1.64 | 2 / 5 |
| convenienceYieldCarry | −0.28 | 1.92 | 1 / 5 |
| storageCashCarry | −1.95 | 3.23 | 1 / 5 |
| **curveRelativeValue** | **−1.76** | 1.71 | **0 / 5** |

The single-split "winner" (curveRelativeValue, +2.19) is **negative in every one
of the five splits** (mean −1.76): the 2021 result was one lucky draw. The
robust performer is **twoTimescale** (revert the fast χ factor + trend-follow the
slow ξ factor), positive in 4 of 5 splits — and even it has a dispersion (1.28)
comparable to its mean, i.e. a wide error bar.

**The honest read** — and the negative findings are the valuable ones:
- A one-to-few-year Sharpe has a standard error of order ±1-2 (visible in the
  dispersions above), so a point estimate is not a verdict; the *distribution* is.
- The cross-sectional **relative-value** edge was a single-draw artifact, not a
  repeatable signal — exactly the trap walk-forward exists to catch.
- **Momentum** being positive in a trending year is partly mechanical, and it was
  positive in only 2 of 5 splits — regime-dependent, not a free lunch.
- **Carry** underperformed because its convenience-yield signal leans on a Kalman
  κ that **railed high on the COVID-distorted 2020 training window** — a model is
  only as good as the regime it is trained on.
- **chi-reversion** didn't trade at all in the single calm 2021 split, yet across
  the volatile rolling splits it was positive in 4 of 5: the signal needs a regime
  with dislocations to act on.
- The two-timescale combination being the most consistent is the one mild
  positive: combining a fast mean-reverting and a slow trend factor diversifies
  the regime dependence of either alone.

**Limitations & caveats (stated plainly).** This is a *single commodity* (WTI) over
a *finite liquid history* dominated by one extraordinary regime (the 2020 COVID
shock); the lognormal models structurally **exclude the April-2020 negative-price
window**; single-curve calibration **fixes the volatilities** (only the curve-
shaping parameters are identified from one snapshot); convergence trades look
flatter and more attractive in-sample than they trade net of cost; and five
overlapping rolling splits are not five independent experiments. The framework's
contribution is the *discipline* — causal signals, out-of-sample walk-forward,
common vol target, costs, and an accounting identity checked by independent
revaluation — not a claim of a money-making strategy.

### Cross-commodity robustness (CL + NG)

One commodity over one window is still anecdote. Extending the WTI walk-forward
to **all 1,866 liquid dates (2018-12 → 2026-05, 10 rolling 252/63 splits)** and
adding **Henry Hub gas** — with its signals **deseasonalized** (monthly factors
fit train-only, then divided out before extracting carry/Kalman signals; the
tradeable returns stay raw) — gives **20 (commodity × split) cells** per strategy:

| strategy | mean OOS Sharpe | dispersion | fraction of cells positive |
|---|---|---|---|
| **timeSeriesMomentum** | **+0.74** | 1.13 | **0.75** |
| carryMomentumCombo | +0.66 | 1.39 | 0.70 |
| chiReversion | +0.51 | 2.25 | 0.60 |
| convenienceYieldCarry | +0.41 | 1.67 | 0.50 |
| twoTimescale | −0.05 | 2.14 | 0.60 |
| curveRelativeValue | −0.31 | 2.00 | 0.50 |
| storageCashCarry | −0.78 | 2.11 | 0.30 |

**What survived, and what didn't:**
- **`twoTimescale` did NOT survive.** It was the §8 walk-forward standout on the
  2020-2021 WTI window (mean +1.77, 4/5) — but over the full 2018-2026 WTI history
  it is **negative** (mean −0.15) and cross-commodity it is **−0.05**. Its v0.54
  "edge" was the COVID-recovery window, not a repeatable signal. Extending the
  window is what exposed it.
- **`timeSeriesMomentum` is the robust performer**: positive mean across both
  commodities and **positive in 75% of cells** — consistent with the documented
  commodity trend premium, the one signal that holds up across commodity *and*
  window. `carryMomentumCombo` (carry + momentum) is the runner-up (0.70 positive).
- **Deseasonalization changed the gas-carry verdict.** Raw NG carry just trades
  the mechanical winter/summer calendar — its out-of-sample Sharpe on a mid-history
  split is **−0.02** (no edge). Deseasonalized (removing the fitted monthly premium
  first), the same signal is **+0.18** on that split and **+0.80** mean across the
  NG walk-forward: only after deseasonalizing does gas carry mean anything.

**Caveat, stated plainly.** This is still **two commodities over ~2018-2026**, a
window dominated by one extraordinary regime (2020). Dispersions of 1-2 Sharpe
units mean even the "robust" momentum result has wide error bars. True multi-regime
robustness needs pre-2020 history (and more commodities); the honest conclusion is
*momentum is the least fragile of these signals here*, not that it is a proven edge.

## 9. Engineering rigor

A few properties make the numbers above trustworthy rather than merely present:

- **Independent-revaluation accounting.** Every strategy's P&L is validated by rebuilding the
  portfolio value from atomic state (`cash + leg marks + hedge·spot`) and asserting it equals
  the reported `stepPnl` — a different computation, not a restatement of the same formula.
  Residuals are at machine epsilon (≤ ~4e-15) across all strategies.
- **Byte-identical guarantees.** Every milestone preserves prior behavior exactly: opt-in
  features (e.g. seasonality) leave the existing adapters byte-for-byte unchanged, and the
  canonical reference tests keep their reference numbers.
- **364 tests, all green**, run in ~6-7 min — kept under a 10-minute ceiling by deliberate
  test-suite-performance work (per-test timing, tiered runners, and grid/path tuning that
  preserves the machine-epsilon identities).
- **Honest validation.** Where a method has a limitation (single-curve κ under-identification),
  it is stated and then *addressed* by a better method, not hidden.

## 10. How to reproduce

These read the real (gitignored) WTI CSVs and print the numbers quoted above:

| script | prints |
|---|---|
| `examples/load_crude_curve.q` | the Jan-2020 forward curve + contango/backwardation read |
| `examples/calibrate_crude_curve.q` | the two-factor fit, RMSE, and convenience-yield-vs-rate read |
| `examples/convenience_yield_series.q` | the 2020-2021 cy series + regime-transition dates |
| `examples/kalman_schwartz_smith.q` | the Kalman-MLE params (κ identified) + filtered χ/ξ factors |
| `examples/commodity_strategy_backtest.q` | the single-split (train 2020 / test 2021) ranked OOS table + correlations |

The walk-forward distribution in §8 comes from `.strategy.commodityBT.walkForward` over
rolling splits (run inside the backtest harness; see `lib/commodityStrategies.q`). Run any
example from the repo root, e.g. `q examples/kalman_schwartz_smith.q`. The synthetic
identifiability round-trip lives in `tests/commodity/test_kalman_schwartz_smith_roundtrip.q`.

---

*The arc: prove the numerics on equity vanillas → build the commodity term-structure models →
validate the spread pricers against each other → ingest real WTI → read the curve's
convenience yield → identify the dynamics the snapshot couldn't → and watch the model narrate
the 2020 oil shock. Each step is checked against something independent.*
