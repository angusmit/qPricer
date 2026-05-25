# qFDM: kdb+/q Finite Difference Pricing Framework

qFDM is a modular kdb+/q pricing framework for equity options using finite-difference methods, with Greeks calculation and scenario risk reporting.

## Current Scope (v0.3)

**Supported:**
- European call and put options
- **American put option** (early exercise via explicit FDM)
- Black-Scholes PDE with explicit finite-difference method
- Flat spot / rate / dividend / volatility inputs
- Greeks: delta, gamma, theta (grid), vega and rho (bump-and-reprice)
- Scenario risk reports (spot, vol, rate bumps)
- Analytical Black-Scholes validation (price and Greeks)
- Put-call parity, grid convergence, input validation

**Not yet supported:**
- Barrier options
- Local volatility / volatility surfaces
- Implicit FDM, Crank-Nicolson
- Portfolio pricing

## Quick Start

```q
\l lib/init.q

/ European call
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);

/ American put — just change exerciseStyle and optionType
americanPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`american;`put;100f;1f;1f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f)];

show .engine.priceOption[trade;marketData;model;config]
show .engine.priceOption[americanPut;marketData;model;config]
```

## Public API

| Function | Returns |
|----------|---------|
| `.engine.priceOption[trade;mkt;model;cfg]` | Price dictionary |
| `.engine.priceOptionWithGrid[trade;mkt;model;cfg]` | Price + full grid |
| `.greeks.calculateGreeks[trade;mkt;model;cfg]` | Greeks table |
| `.risk.generateScenarioReport[trade;mkt;model;cfg]` | Scenario risk table |
| `.validation.validateEuropeanOption[trade;mkt;model;cfg]` | FDM vs BS price |
| `.validation.validateGreeks[trade;mkt;model;cfg]` | FDM vs BS Greeks |
| `.validation.checkPutCallParity[call;put;mkt;model;cfg]` | Parity check |
| `.validation.runGridConvergenceTest[trade;mkt;model;cfgList]` | Convergence table |

## American Put (v0.3)

American puts can be exercised at any time before expiry. The explicit FDM handles this by applying, at each backward time step:

```
optionValue = max(continuationValue, intrinsicValue)
```

where `intrinsicValue = max(K - S, 0)` for a put.

### Usage

To price an American put, set `exerciseStyle:\`american` and `optionType:\`put`:

```q
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`put;100f;1f;1f);
show .engine.priceOption[trade;marketData;model;config]
```

### Early Exercise Premium

The American put price is always >= the European put price. The difference is the early exercise premium, which is positive when interest rates are positive (exercising early captures the time value of the strike proceeds).

For ATM options (S=100, K=100, T=1y, r=5%, vol=20%):

| Metric | Value |
|--------|-------|
| European put | ~5.578 |
| American put | ~5.8-6.1 (grid dependent) |
| Early exercise premium | > 0 |

### Limitations

- Only American **put** has meaningful early exercise premium with zero dividends
- American call on non-dividend stock = European call (no premium)
- Uses explicit FDM only (Crank-Nicolson would be more accurate, planned for future)

## Scenario Risk Reports

9 scenarios: base, spot ±1%/±5%, vol ±1pp, rate ±25bp. See v1.2 documentation.

## Greeks

Delta, gamma, theta from FDM grid. Vega, rho by bump-and-reprice. See v1.1 documentation.

## Validation

### European Price (200 spot, 2000 time, [0,300])

| Option | FDM | BS | Error |
|--------|-----|----|-------|
| Call | 10.45496 | 10.45058 | 0.00438 |
| Put | 5.57784 | 5.57352 | 0.00432 |

## Architecture

```
lib/
  init.q          Silent loader
  utilities.q     Validation, interpolation
  config.q        FDM configuration
  product.q       Trade definition (european + american)
  market.q        Market data + bumping
  model.q         Black-Scholes model
  grid.q          Grid construction
  payoff.q        Terminal payoff
  boundary.q      Boundary conditions
  solver.q        Explicit FDM solver + early exercise
  engine.q        Pricing engine
  greeks.q        Greeks (grid + bump)
  validation.q    BS closed form + validation
  risk.q          Scenario risk reports
```

## Running

```
q examples/smoke_test_european_call.q       # European call
q examples/price_american_put.q             # American put vs European
q examples/calculate_greeks.q               # Greeks
q examples/generate_scenario_report.q       # Scenario risk
q tests/run_all_tests.q                     # full suite (9 tests)
```

## Future Extensions

- Barrier option boundary conditions
- Implicit FDM / Crank-Nicolson (better for American options)
- Local volatility sigma(S,t)
- Volatility surface input
- American call with dividends
- Portfolio-level pricing
- HDB market data integration
