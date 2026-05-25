# qFDM: kdb+/q Finite Difference Pricing Framework

qFDM is a modular kdb+/q pricing framework for equity options using finite-difference methods.

## Current Scope (v0.6)

**Supported:**
- European call and put (explicit + Crank-Nicolson)
- American put (explicit only)
- Up-and-out call, down-and-out put (explicit only)
- Greeks, scenario risk, exercise boundary extraction
- Black-Scholes analytical validation

**Not yet supported:**
- CN for American/barrier options
- Implicit method, local volatility, portfolio pricing

## Quick Start

```q
\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];

/ Explicit
explicitCfg:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f)];

/ Crank-Nicolson (fewer time steps needed)
cnCfg:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `crankNicolson;200;500;300f)];

show .engine.priceOption[trade;marketData;model;explicitCfg]
show .engine.priceOption[trade;marketData;model;cnCfg]
```

## Crank-Nicolson Solver (v0.6)

CN averages explicit and implicit time stepping, requiring a tridiagonal solve at each step. It is unconditionally stable, so fewer time steps are needed compared to the explicit method.

### Supported in v0.6
- European call and put under Black-Scholes with flat market data

### Not supported in v0.6
- American options, barrier options, local volatility

### Usage
Set `method:\`crankNicolson` in config. Everything else is the same API.

### Comparison
| Method | Spot Steps | Time Steps | Call Price | BS Price | Error |
|--------|-----------|-----------|-----------|---------|-------|
| Explicit | 200 | 2000 | 10.455 | 10.451 | 0.004 |
| CN | 200 | 500 | ~10.451 | 10.451 | ~0.001 |

CN achieves comparable or better accuracy with 4x fewer time steps.

## Public API

| Function | Returns |
|----------|---------|
| `.engine.priceOption[trade;mkt;model;cfg]` | Price dictionary |
| `.engine.priceOptionWithGrid[trade;mkt;model;cfg]` | Price + full grid |
| `.greeks.calculateGreeks[trade;mkt;model;cfg]` | Greeks table |
| `.risk.generateScenarioReport[trade;mkt;model;cfg]` | Scenario risk table |
| `.american.extractEarlyExerciseBoundary[trade;mkt;model;cfg]` | Exercise boundary |
| `.american.analyzeAmericanPut[trade;mkt;model;cfg]` | Price + premium + boundary |
| `.validation.validateEuropeanOption[trade;mkt;model;cfg]` | FDM vs BS price |
| `.validation.validateGreeks[trade;mkt;model;cfg]` | FDM vs BS Greeks |

## Architecture

```
lib/
  init.q          Silent loader
  utilities.q     Validation, interpolation
  config.q        FDM configuration (explicit + CN)
  product.q       Trade definition (vanilla + barrier)
  market.q        Market data + bumping
  model.q         Black-Scholes model
  grid.q          Grid construction
  payoff.q        Terminal payoff
  boundary.q      European boundaries + barrier
  solver.q        Explicit + Crank-Nicolson + Thomas algorithm
  engine.q        Pricing engine (method routing)
  greeks.q        Greeks
  validation.q    BS closed form + validation
  risk.q          Scenario risk
  american.q      American put analysis
```

## Running

```
q examples/smoke_test_european_call.q       # European call
q examples/compare_explicit_crank_nicolson.q # Explicit vs CN
q examples/price_american_put.q             # American put
q examples/price_barrier_options.q          # Barrier options
q examples/calculate_greeks.q               # Greeks
q examples/generate_scenario_report.q       # Scenario risk
q tests/run_all_tests.q                     # full suite (16 tests)
```

## Version History

| Version | Features |
|---------|----------|
| v0.1 | European call/put, explicit FDM, BS validation |
| v0.2 | Greeks, scenario risk, input validation |
| v0.3 | American put, early exercise premium |
| v0.4 | Early exercise boundary extraction |
| v0.5 | Up-and-out call, down-and-out put |
| v0.6 | Crank-Nicolson for European vanilla |

## Roadmap

| Version | Planned |
|---------|---------|
| v0.7 | CN for American options, or local volatility |
| v0.8 | Knock-in barriers, rebates |
| v0.9 | Portfolio-level pricing |
