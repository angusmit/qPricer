# qFDM: kdb+/q Finite Difference Pricing Framework

qFDM is a modular kdb+/q pricing framework for European equity options using finite-difference methods, with Greeks calculation and scenario risk reporting.

## Current Scope (v0.2)

**Supported:**
- European call and put options
- Black-Scholes PDE with explicit finite-difference method
- Flat spot / rate / dividend / volatility inputs
- Greeks: delta, gamma, theta (grid), vega and rho (bump-and-reprice)
- Scenario risk reports (spot, vol, rate bumps)
- Analytical Black-Scholes validation (price and Greeks)
- Put-call parity, grid convergence, input validation

**Not yet supported:**
- American exercise, barrier options
- Local volatility / volatility surfaces
- Implicit FDM, Crank-Nicolson
- Portfolio pricing

## Quick Start

```q
\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1000000f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f)];

show .engine.priceOption[trade;marketData;model;config]
show .greeks.calculateGreeks[trade;marketData;model;config]
show .risk.generateScenarioReport[trade;marketData;model;config]
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

## Scenario Risk Reports

qFDM v0.2 generates risk reports by bumping market data and repricing through the engine.

### Scenarios

| Scenario | Bump Type | Bump Size |
|----------|-----------|-----------|
| `base` | none | — |
| `spotUp1Pct` | spot relative | +1% |
| `spotDown1Pct` | spot relative | -1% |
| `spotUp5Pct` | spot relative | +5% |
| `spotDown5Pct` | spot relative | -5% |
| `volatilityUp1Point` | vol absolute | +1pp |
| `volatilityDown1Point` | vol absolute | -1pp |
| `rateUp25Bp` | rate absolute | +25bp |
| `rateDown25Bp` | rate absolute | -25bp |

### Output columns

`tradeId`, `underlying`, `optionType`, `scenario`, `spot`, `riskFreeRate`, `dividendYield`, `volatility`, `unitPrice`, `unitPnL`, `notionalPrice`, `notionalPnL`

### Conventions

- `unitPnL` = scenario unit price - base unit price
- `notionalPnL` = unitPnL * notional
- Spot bumps are relative (0.01 = +1%)
- Vol bumps are absolute (0.01 = +1 vol point)
- Rate bumps are absolute (0.0025 = +25bp)

## Greeks

### Conventions

| Greek | Convention | Example |
|-------|-----------|---------|
| Delta | Per 1 unit spot move | delta=0.637: +$1 spot -> +$0.637 price |
| Gamma | Per 1 unit spot squared | gamma=0.019 |
| Theta | Annual (dV/dt) | theta=-6.41: loses $6.41/year |
| Vega | Per 1.00 absolute vol | vega=37.5: +1pp vol -> +$0.375 price |
| Rho | Per 1.00 absolute rate | rho=53.2: +1bp rate -> +$0.00532 price |

### Reference (S=100, K=100, T=1y, r=5%, q=0%, vol=20%)

| Greek | Call | Put |
|-------|------|-----|
| Delta | 0.6368 | -0.3632 |
| Gamma | 0.0188 | 0.0188 |
| Theta | -6.414 | -1.658 |
| Vega | 37.524 | 37.524 |
| Rho | 53.233 | -41.891 |

## Validation

### Price (200 spot steps, 2000 time steps, [0,300])

| Option | FDM | BS | Error |
|--------|-----|----|-------|
| Call | 10.45496 | 10.45058 | 0.00438 |
| Put | 5.57784 | 5.57352 | 0.00432 |

### Grid Convergence

| Spot Steps | Time Steps | Error |
|-----------|-----------|-------|
| 50 | 250 | 0.06717 |
| 200 | 2000 | 0.00438 |

## Architecture

```
lib/
  init.q          Silent loader
  utilities.q     Validation, interpolation
  config.q        FDM configuration
  product.q       Trade definition
  market.q        Market data + bumping
  model.q         Black-Scholes model
  grid.q          Grid construction
  payoff.q        Terminal payoff
  boundary.q      Boundary conditions
  solver.q        Explicit FDM solver
  engine.q        Pricing engine
  greeks.q        Greeks (grid + bump)
  validation.q    BS closed form + validation
  risk.q          Scenario risk reports
```

## Running

```
q examples/smoke_test_european_call.q       # price + validate
q examples/calculate_greeks.q               # Greeks + validate
q examples/generate_scenario_report.q       # price + Greeks + risk
q tests/run_all_tests.q                     # full suite (8 tests)
```

## Future Extensions

- American option early exercise
- Barrier option boundary conditions
- Implicit FDM / Crank-Nicolson
- Local volatility sigma(S,t)
- Volatility surface input
- Portfolio-level pricing
- HDB market data integration
- IPC/API access
