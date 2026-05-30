# qFDM

A kdb+/q finite-difference pricing and risk framework for equity options.

kdb+/q finite-difference pricing and risk framework for equity options, including explicit FDM, Crank-Nicolson, Greeks, scenario risk, American puts, barriers, local volatility and portfolio analytics.

> **Commodity narrative:** the equity FDM is the correctness-validation case; the project's target is commodity. See **[docs/COMMODITY_DESK.md](docs/COMMODITY_DESK.md)** for the end-to-end story (WTI curve → two-factor calibration → Kalman state-space estimation → the 2020-2021 convenience-yield regimes) told with real numbers.

## Why This Project Exists

qFDM is a kdb+/q pricing and risk analytics framework for equity options. It implements finite-difference pricing under Black-Scholes, American early exercise, simple knock-out barriers, Crank-Nicolson for European vanilla options, local volatility experiments, Greeks, scenario risk, and portfolio-level batch analytics.

The goal is to demonstrate how numerical pricing, validation, and desk-style risk analytics can be structured in q as a modular library rather than a single pricing script.

## Feature Summary

| Area | Supported |
|------|-----------|
| Products | European call/put, American put, up-and-out call, down-and-out put |
| Methods | Explicit FDM, Crank-Nicolson for European vanilla |
| Models | Black-Scholes, local volatility for European vanilla explicit FDM |
| Risk | Delta, Gamma, Theta, Vega, Rho, scenario risk |
| Portfolio | Batch pricing, portfolio Greeks, portfolio scenario risk |
| Validation | Black-Scholes closed form, put-call parity, grid convergence, Greek validation |
| Tests | 24 passing tests |

### Not Yet Supported

- Dupire calibration
- Real implied volatility surface ingestion
- Stochastic volatility
- Jump diffusion
- XVA
- Multi-asset products
- Portfolio-level multi-symbol market data
- HDB persistence
- IPC/server deployment
- Production trading controls

## Architecture Overview

```
Market Data
    |
    v
Product / Trade Definition
    |
    v
Model
    |
    v
Grid / Payoff / Boundary
    |
    v
Solver
    |
    v
Engine
    |
    v
Greeks / Risk / American / Portfolio
    |
    v
Validation / Examples / Tests
```

- `product.q` defines and validates trades.
- `market.q` provides spot, rates, dividends, volatility and local volatility.
- `model.q` defines Black-Scholes and local-vol model assumptions.
- `grid.q` builds spot/time grids.
- `payoff.q` defines payoff and intrinsic value.
- `boundary.q` applies European and barrier boundary conditions.
- `solver.q` contains explicit and Crank-Nicolson solvers.
- `engine.q` is the clean public pricing API.
- `greeks.q` computes sensitivities.
- `risk.q` generates scenario reports.
- `american.q` extracts early exercise boundary.
- `portfolio.q` performs table-based batch analytics.
- `validation.q` contains analytical checks and regression validation.

## Folder Structure

```
q-fdm-option-pricer/
├── lib/
│   ├── init.q
│   ├── utilities.q
│   ├── config.q
│   ├── product.q
│   ├── market.q
│   ├── model.q
│   ├── grid.q
│   ├── payoff.q
│   ├── boundary.q
│   ├── solver.q
│   ├── engine.q
│   ├── greeks.q
│   ├── validation.q
│   ├── risk.q
│   ├── american.q
│   └── portfolio.q
├── examples/
├── tests/
└── README.md
```

## Quick Start

```q
\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);

model:.model.createBlackScholesModel[];

config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

.engine.priceOption[trade;marketData;model;config]
```

Expected output:

```
tradeId      | 1
underlying   | `AAPL
optionType   | `call
unitPrice    | 10.45496
notionalPrice| 10.45496
method       | `explicit
modelName    | `blackScholes
```

## Single-Trade Pricing Example

```q
priceResult:.engine.priceOption[trade;marketData;model;config]
```

Returns a dictionary with tradeId, underlying, optionType, unitPrice, notionalPrice, method, and modelName.

## Greeks Example

```q
show .greeks.calculateGreeks[trade;marketData;model;config]
```

```
tradeId underlying optionType delta     gamma      theta     vega     rho
1       AAPL       call       0.6177    0.01926    -6.412    37.51    53.21
```

Delta and gamma are from the FDM grid via central differences. Vega and rho are by bump-and-reprice.

## Scenario Risk Example

```q
show .risk.generateScenarioReport[trade;marketData;model;config]
```

Returns 9 rows: base, spotUp1Pct, spotDown1Pct, spotUp5Pct, spotDown5Pct, volatilityUp1Point, volatilityDown1Point, rateUp25Bp, rateDown25Bp. Each row includes unitPrice, unitPnL, notionalPrice, and notionalPnL.

## American Put Example

```q
americanPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`american;`put;100f;1f;1f);

show .engine.priceOption[americanPut;marketData;model;config]
/ unitPrice: 6.093274 (vs European put 5.577839, early exercise premium 0.5154354)

show .american.extractEarlyExerciseBoundary[americanPut;marketData;model;config]
/ Table with timePoint, remainingTime, exerciseBoundary per time step
```

## Barrier Option Example

```q
barrierCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    3;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;0f);

show .engine.priceOption[barrierCall;marketData;model;config]
/ unitPrice: 3.330 (vs vanilla call 10.455, discounted by barrier)
```

Barrier condition is applied directly on the FDM grid: option value set to zero at and beyond the barrier level at each time step.

## Local Volatility Example

```q
/ Downside-only local vol uplift: vol >= 0.2 everywhere
localVolFn:{[spotValue;timePoint] 0.2 + (0f | 0.001 * (100f - spotValue))};
localVolMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;localVolFn];
lvModel:.model.createLocalVolatilityModel[];

show .engine.priceOption[trade;localVolMkt;lvModel;config]
/ unitPrice: 10.567 (higher than flat 10.455 due to downside vol uplift)
```

A constant local volatility function reproduces the flat Black-Scholes explicit FDM price exactly (diff = 0). Non-flat local vol functions change prices, confirming the solver correctly uses sigma(S,t) at each grid point.

## Portfolio Pricing Example

```q
tradeTable:([]
    tradeId:1 2 3 4;
    underlying:`AAPL`AAPL`AAPL`AAPL;
    productType:`equityOption`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`american`european;
    optionType:`call`put`put`call;
    strike:100 100 100 100f;
    expiry:1 1 1 1f;
    notional:1000000 1000000 1000000 1000000f;
    barrierType:`none`none`none`upAndOut;
    barrierLevel:0N 0N 0N 130f;
    rebate:0 0 0 0f);

show .portfolio.priceTrades[tradeTable;marketData;model;config]
/ 4 rows: one per trade with unitPrice, notionalPrice, status

show .portfolio.calculatePortfolioGreeks[tradeTable;marketData;model;config]
/ European vanilla: Greeks populated. American/barrier: status `UNSUPPORTED.

show .portfolio.generatePortfolioScenarioReport[tradeTable;marketData;model;config]
/ 36 rows (4 trades x 9 scenarios), with unitPnL and notionalPnL
```

Each trade is priced independently with per-trade error handling. A failed trade returns status `ERROR without crashing the batch.

## Validation and Test Suite

```
qFDM v0.8.1 Test Suite

Results: 24 passed, 0 failed
All tests passed.
```

### Test Categories

- European call validation
- European put validation
- Put-call parity
- Grid convergence
- Input validation (23 edge cases)
- Greeks call/put validation
- Scenario risk
- American put premium
- Early exercise boundary
- Barrier option pricing (up-and-out call, down-and-out put)
- Barrier validation (7 invalid combinations)
- Crank-Nicolson call/put/vs-explicit
- Local volatility flat equivalence call/put
- Local volatility validation (6 unsupported combinations)
- Local volatility skew sanity
- Local volatility time dependence
- Portfolio pricing
- Portfolio Greeks
- Portfolio scenario risk

### Sample Test Output

European call explicit:
- FDM = 10.45496
- Black-Scholes = 10.45058
- error = 0.004379982

European put explicit:
- FDM = 5.577839
- Black-Scholes = 5.573518
- error = 0.004320529

Crank-Nicolson call:
- CN = 10.45443
- Black-Scholes = 10.45058
- error = 0.003857219

American put:
- American = 6.093274
- European = 5.577839
- early exercise premium = 0.5154354

Barrier examples:
- up-and-out call = 3.32987 vs vanilla call = 10.44837
- down-and-out put = 4.00776 vs vanilla put = 5.571288

Local volatility:
- flat local vol equivalence diff = 0
- downside-only local vol uplift increases call and put prices

Portfolio:
- 4-trade portfolio pricing passes
- portfolio scenario risk returns 36 rows
- portfolio Greeks returns supported and unsupported rows correctly

## How This Maps to Trading Desk Workflows

qFDM is not a full production trading system. It is a simplified educational/research framework showing the structure of desk analytics:

- Price a trade - `.engine.priceOption`
- Validate model output - `.validation.validateEuropeanOption`
- Calculate Greeks - `.greeks.calculateGreeks`
- Run scenario shocks - `.risk.generateScenarioReport`
- Compare American vs European early exercise value - `.american.analyzeAmericanPut`
- Monitor barrier behaviour - barrier trade pricing under different scenarios
- Batch-price a portfolio table - `.portfolio.priceTrades`
- Generate portfolio scenario risk - `.portfolio.generatePortfolioScenarioReport`

In a real desk system, market data, curves, vol surfaces, trade capture, risk aggregation, limits, and production controls would be much more complex. qFDM focuses on the numerical pricing/risk core and q-style table analytics.

## Module Responsibility

| Module | Responsibility |
|--------|----------------|
| `config.q` | Finite-difference configuration and validation |
| `product.q` | Trade dictionaries and product validation |
| `market.q` | Spot/rate/dividend/volatility/local volatility market data |
| `model.q` | Model selection: Black-Scholes and local volatility |
| `grid.q` | Spot/time grid construction |
| `payoff.q` | Payoff and intrinsic value |
| `boundary.q` | European and barrier boundary conditions |
| `solver.q` | Explicit FDM and Crank-Nicolson solvers |
| `engine.q` | Public pricing API |
| `greeks.q` | Delta/Gamma/Theta/Vega/Rho |
| `risk.q` | Single-trade scenario risk |
| `american.q` | American early exercise analysis |
| `portfolio.q` | Portfolio pricing and batch risk |
| `validation.q` | Analytical validation and regression checks |

## Public API

| Function | Returns |
|----------|---------|
| `.engine.priceOption[trade;mkt;model;cfg]` | Price dictionary |
| `.engine.priceOptionWithGrid[trade;mkt;model;cfg]` | Price + full grid |
| `.greeks.calculateGreeks[trade;mkt;model;cfg]` | Greeks table |
| `.risk.generateScenarioReport[trade;mkt;model;cfg]` | Scenario risk table |
| `.american.extractEarlyExerciseBoundary[trade;mkt;model;cfg]` | Exercise boundary table |
| `.american.analyzeAmericanPut[trade;mkt;model;cfg]` | Price + premium + boundary |
| `.portfolio.priceTrades[tbl;mkt;model;cfg]` | Portfolio prices |
| `.portfolio.calculatePortfolioGreeks[tbl;mkt;model;cfg]` | Portfolio Greeks |
| `.portfolio.generatePortfolioScenarioReport[tbl;mkt;model;cfg]` | Portfolio scenarios |
| `.portfolio.summarizePortfolioRisk[scenarioTbl]` | Scenario summary |
| `.validation.validateEuropeanOption[trade;mkt;model;cfg]` | FDM vs BS |
| `.validation.validateGreeks[trade;mkt;model;cfg]` | FDM vs BS Greeks |
| `.validation.checkPutCallParity[call;put;mkt;model;cfg]` | Parity check |
| `.validation.runGridConvergenceTest[trade;mkt;model;cfgList]` | Convergence table |

## Running

```
cd q-fdm-option-pricer

q tests/run_all_tests.q
q examples/smoke_test_european_call.q
q examples/calculate_greeks.q
q examples/generate_scenario_report.q
q examples/price_american_put.q
q examples/analyze_american_put.q
q examples/price_barrier_options.q
q examples/compare_explicit_crank_nicolson.q
q examples/price_local_volatility.q
q examples/price_local_volatility_skew.q
q examples/price_portfolio.q
```

## Version History

| Version | Features |
|---------|----------|
| v0.1 | European call/put explicit FDM, Black-Scholes validation, put-call parity, grid convergence |
| v0.2 | Greeks, scenario risk, input validation |
| v0.3 | American put pricing and early exercise premium |
| v0.4 | Early exercise boundary extraction |
| v0.5 | Up-and-out call and down-and-out put |
| v0.6 | Crank-Nicolson solver for European vanilla options |
| v0.7 | Local volatility for European vanilla explicit FDM |
| v0.7.1 | Local volatility sanity checks and downside-only uplift example |
| v0.8 | Portfolio pricing and batch risk |
| v0.8.1 | Documentation and project polish |

## Future Roadmap

- v0.9: market data ingestion / multi-symbol market data
- v1.0: cleaned release with docs, examples, benchmarks
- Later: Dupire local volatility calibration
- Later: HDB storage
- Later: IPC API
- Later: performance benchmarking
- Later: CI setup

## Limitations

- This is not a production trading system.
- The local volatility examples are toy functions, not calibrated Dupire surfaces.
- Barrier pricing uses simple grid knock-out conditions and no rebate support.
- Crank-Nicolson currently supports only European vanilla options.
- Portfolio module currently assumes one shared market data object.
- No multi-underlying market data routing yet.
- No HDB/persistence yet.
- No parallel pricing yet.
- No CI setup yet.
- Explicit FDM has stability constraints.
- Greeks are approximate and method-dependent.

## Disclaimer

This project is for educational and demonstration purposes. It is not production trading software and should not be used for actual financial decisions.