# qFDM: kdb+/q Finite Difference Pricing Framework

qFDM is a modular kdb+/q pricing framework for equity options using finite-difference methods, with Greeks, scenario risk, and American option support.

## Current Scope (v0.4)

**Supported:**
- European call and put
- American put (early exercise via explicit FDM)
- American put early exercise boundary extraction
- Black-Scholes PDE with explicit finite-difference method
- Greeks: delta, gamma, theta (grid), vega and rho (bump-and-reprice)
- Scenario risk reports (spot, vol, rate bumps)
- Analytical Black-Scholes validation (price and Greeks)
- Put-call parity, grid convergence, input validation

**Not yet supported:**
- American call with dividends
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

/ American put - just change exerciseStyle and optionType
americanPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`american;`put;100f;1f;1f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f)];

show .engine.priceOption[trade;marketData;model;config]
show .engine.priceOption[americanPut;marketData;model;config]
show .american.extractEarlyExerciseBoundary[americanPut;marketData;model;config]
```

## Public API

| Function | Returns |
|----------|---------|
| `.engine.priceOption[trade;mkt;model;cfg]` | Price dictionary |
| `.engine.priceOptionWithGrid[trade;mkt;model;cfg]` | Price + full grid |
| `.greeks.calculateGreeks[trade;mkt;model;cfg]` | Greeks table |
| `.risk.generateScenarioReport[trade;mkt;model;cfg]` | Scenario risk table |
| `.american.extractEarlyExerciseBoundary[trade;mkt;model;cfg]` | Exercise boundary table |
| `.american.analyzeAmericanPut[trade;mkt;model;cfg]` | Price + premium + boundary |
| `.validation.validateEuropeanOption[trade;mkt;model;cfg]` | FDM vs BS price |
| `.validation.validateGreeks[trade;mkt;model;cfg]` | FDM vs BS Greeks |
| `.validation.checkPutCallParity[call;put;mkt;model;cfg]` | Parity check |
| `.validation.runGridConvergenceTest[trade;mkt;model;cfgList]` | Convergence table |

## American Put (v0.3-v0.4)

### Early Exercise

American puts can be exercised at any time. The explicit FDM applies at each backward time step:

```
optionValue = max(continuationValue, intrinsicValue)
```

### Early Exercise Premium

The American put price always exceeds the European put price:

| Metric | Value |
|--------|-------|
| American put | 6.093 |
| European put | 5.578 |
| Early exercise premium | 0.515 |

(S=100, K=100, T=1y, r=5%, q=0%, vol=20%, 200 spot steps, 2000 time steps)

### Early Exercise Boundary (v0.4)

The boundary is the highest spot where immediate exercise is optimal at each time. Extracted from the FDM grid by comparing option value with intrinsic value.

```q
boundaryTable:.american.extractEarlyExerciseBoundary[trade;marketData;model;config]
```

Returns: `tradeId`, `underlying`, `optionType`, `timePoint`, `remainingTime`, `exerciseBoundary`, `hasExerciseRegion`

The boundary typically starts below the strike near valuation and approaches the strike at expiry.

## Scenario Risk Reports

9 scenarios: base, spot +/-1%/+/-5%, vol +/-1pp, rate +/-25bp.

## Greeks

Delta, gamma, theta from FDM grid. Vega, rho by bump-and-reprice.

| Greek | Convention |
|-------|-----------|
| Delta | Per 1 unit spot move |
| Gamma | Per 1 unit spot squared |
| Theta | Annual (dV/dt) |
| Vega | Per 1.00 absolute vol |
| Rho | Per 1.00 absolute rate |

## Validation

### European Price

| Option | FDM | BS | Error |
|--------|-----|----|-------|
| Call | 10.455 | 10.451 | 0.004 |
| Put | 5.578 | 5.574 | 0.004 |

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
  american.q      American put analysis + exercise boundary
```

## Running

```
q examples/smoke_test_european_call.q       # European call
q examples/price_american_put.q             # American put vs European
q examples/analyze_american_put.q           # Exercise boundary analysis
q examples/calculate_greeks.q               # Greeks
q examples/generate_scenario_report.q       # Scenario risk
q tests/run_all_tests.q                     # full suite (10 tests)
```

## Version History

| Version | Features |
|---------|----------|
| v0.1 | European call/put, explicit FDM, BS validation, put-call parity, grid convergence |
| v0.2 | Greeks, scenario risk report, input validation |
| v0.3 | American put pricing, early exercise premium |
| v0.4 | Early exercise boundary extraction |

## Roadmap

| Version | Planned |
|---------|---------|
| v0.5 | Crank-Nicolson or barrier options |
| v0.6 | Local volatility / volatility surfaces |
| v0.7 | Portfolio-level pricing |

## Future Extensions

- Barrier option boundary conditions
- Implicit FDM / Crank-Nicolson
- Local volatility sigma(S,t)
- Volatility surface input
- American call with dividends
- Portfolio-level pricing
- HDB market data integration
- IPC/API access
