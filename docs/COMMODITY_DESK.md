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

## 8. Engineering rigor

A few properties make the numbers above trustworthy rather than merely present:

- **Independent-revaluation accounting.** Every strategy's P&L is validated by rebuilding the
  portfolio value from atomic state (`cash + leg marks + hedge·spot`) and asserting it equals
  the reported `stepPnl` — a different computation, not a restatement of the same formula.
  Residuals are at machine epsilon (≤ ~4e-15) across all strategies.
- **Byte-identical guarantees.** Every milestone preserves prior behavior exactly: opt-in
  features (e.g. seasonality) leave the existing adapters byte-for-byte unchanged, and the
  canonical reference tests keep their reference numbers.
- **345 tests, all green**, run in ~431 s — kept under a 10-minute ceiling by deliberate
  test-suite-performance work (per-test timing, tiered runners, and grid/path tuning that
  preserves the machine-epsilon identities).
- **Honest validation.** Where a method has a limitation (single-curve κ under-identification),
  it is stated and then *addressed* by a better method, not hidden.

## 9. How to reproduce

All four read the real (gitignored) WTI CSVs and print the numbers quoted above:

| script | prints |
|---|---|
| `examples/load_crude_curve.q` | the Jan-2020 forward curve + contango/backwardation read |
| `examples/calibrate_crude_curve.q` | the two-factor fit, RMSE, and convenience-yield-vs-rate read |
| `examples/convenience_yield_series.q` | the 2020-2021 cy series + regime-transition dates |
| `examples/kalman_schwartz_smith.q` | the Kalman-MLE params (κ identified) + filtered χ/ξ factors |

Run any of them from the repo root, e.g. `q examples/kalman_schwartz_smith.q`. The synthetic
identifiability round-trip lives in `tests/commodity/test_kalman_schwartz_smith_roundtrip.q`.

---

*The arc: prove the numerics on equity vanillas → build the commodity term-structure models →
validate the spread pricers against each other → ingest real WTI → read the curve's
convenience yield → identify the dynamics the snapshot couldn't → and watch the model narrate
the 2020 oil shock. Each step is checked against something independent.*
