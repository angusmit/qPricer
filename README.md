# qFDM: kdb+/q Finite Difference Pricing Framework

qFDM is a modular kdb+/q finite-difference pricing framework for European equity options, structured like a simplified front-office analytics library. It separates product definition, market data, model assumptions, grid construction, solver logic, pricing engine, Greeks, and validation into independent modules.

## Current v1 Scope

**Supported:**
- European call and put options
- Black-Scholes PDE
- Explicit finite-difference method
- Flat spot / rate / dividend / volatility inputs
- Linear interpolation from grid
- Greeks: delta, gamma, theta (grid), vega and rho (bump-and-reprice)
- Black-Scholes closed-form validation
- Put-call parity check
- Grid convergence testing
- Input validation with clear error messages

**Not yet supported:**
- American exercise
- Barrier options
- Local volatility / volatility surfaces
- Implicit finite difference
- Crank-Nicolson
- Portfolio pricing
- Scenario risk reports

## Quick Start

```q
\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);

model:.model.createBlackScholesModel[];

config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f)];

show .engine.priceOption[trade;marketData;model;config]
```

## Public API

| Function | Returns |
|----------|---------|
| `.engine.priceOption[trade;mkt;model;cfg]` | Price dictionary |
| `.engine.priceOptionWithGrid[trade;mkt;model;cfg]` | Price + full grid |
| `.greeks.calculateGreeks[trade;mkt;model;cfg]` | Greeks table |
| `.validation.validateEuropeanOption[trade;mkt;model;cfg]` | FDM vs BS comparison |
| `.validation.checkPutCallParity[call;put;mkt;model;cfg]` | Parity check |
| `.validation.runGridConvergenceTest[trade;mkt;model;cfgList]` | Convergence table |

## Architecture

```
lib/
  init.q          Silent loader, sets .qfdm.loaded
  utilities.q     Validation, interpolation, error helpers
  config.q        FDM configuration with defaults
  product.q       Trade definition and validation
  market.q        Market data: creation, access, bumping
  model.q         Model definition (Black-Scholes)
  grid.q          Spot and time grid construction
  payoff.q        Terminal payoff
  boundary.q      European boundary conditions
  solver.q        Explicit FDM solver (vectorised)
  engine.q        Pricing engine
  greeks.q        Greeks calculation
  validation.q    BS closed form and validation
```

### Naming Conventions

- `.namespace.publicFunction` — intended public API
- `.namespace.__internalHelper` — double underscore = internal

## Validation Results

Test parameters: S=100, K=100, T=1y, r=5%, q=0%, vol=20%, spotSteps=200, timeSteps=2000, spotRange=[0,300]

### European Call

| Metric | Value |
|--------|-------|
| FDM price | 10.45496 |
| Black-Scholes price | 10.45058 |
| Absolute error | 0.00438 |
| Relative error | 0.042% |

### European Put

| Metric | Value |
|--------|-------|
| FDM price | 5.57784 |
| Black-Scholes price | 5.57352 |
| Absolute error | 0.00432 |
| Relative error | 0.078% |

### Put-Call Parity

| Metric | Value |
|--------|-------|
| Actual (C - P) | 4.87712 |
| Theoretical | 4.87706 |
| Absolute error | 0.00006 |

### Grid Convergence

| Spot Steps | Time Steps | Absolute Error |
|-----------|-----------|----------------|
| 50 | 250 | 0.06717 |
| 200 | 2000 | 0.00438 |

As the grid is refined, the explicit FDM price converges toward the Black-Scholes closed-form price.

## Grid Convention

- **Rows** = spot price levels (ascending, 0 to maximumSpot)
- **Columns** = time steps (ascending, t=0 to t=expiry)
- **Column 0** = valuation date
- **Last column** = expiry

The solver steps backward internally, then reverses column order before returning.

## Running Tests

Smoke test:
```
q examples/smoke_test_european_call.q
```

Full test suite:
```
q tests/run_all_tests.q
```

Expected output:
```
--- Running: tests/test_european_call.q ---
PASS test_european_call: FDM=10.45496, BS=10.45058, error=0.004379982

--- Running: tests/test_european_put.q ---
PASS test_european_put: FDM=5.577839, BS=5.573518, error=0.004320529

--- Running: tests/test_put_call_parity.q ---
PASS test_put_call_parity: actual=4.877117, theoretical=4.877058, error=5.945281e-005

--- Running: tests/test_grid_convergence.q ---
PASS test_grid_convergence: coarse error=0.06717266, fine error=0.004379982

--- Running: tests/test_input_validation.q ---
PASS test_input_validation

Results: 5 passed, 0 failed
All tests passed.
```

## Future Extensions

- American option early exercise
- Barrier option boundary conditions
- Implicit finite difference (tridiagonal solver)
- Crank-Nicolson (second-order in time)
- Local volatility sigma(S,t)
- Volatility surface input
- Risk scenario reports
- Portfolio-level pricing
- HDB market data integration
- IPC/API access from another q process
