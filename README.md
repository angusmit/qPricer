# qFDM: kdb+/q Finite Difference Pricing Framework

qFDM is a modular kdb+/q pricing framework for equity options using finite-difference methods.

## Current Scope (v0.8)

**Supported:**
- European call and put (explicit + Crank-Nicolson)
- American put (explicit only)
- Up-and-out call, down-and-out put (explicit only)
- Local volatility sigma(S,t) for European vanilla (explicit only)
- Portfolio batch pricing, Greeks, and scenario risk
- Greeks, scenario risk, exercise boundary extraction
- Black-Scholes analytical validation

**Not yet supported:**
- Multi-underlying portfolios, parallel pricing
- CN/American/barrier with local vol
- Dupire calibration, implied vol surface
- HDB, IPC, persistence

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

## Portfolio Pricing (v0.8)

A portfolio is a q table of trades. Each row is priced through the single-trade engine with per-trade error handling.

### Trade Table

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
    barrierLevel:0Nf 0Nf 0Nf 130f;
    rebate:0 0 0 0f);
```

### Batch Pricing

```q
portfolioPrices:.portfolio.priceTrades[tradeTable;marketData;model;config]
```

Returns one row per trade with unitPrice, notionalPrice, status, statusMessage. Failed trades return status:\`ERROR without crashing the batch.

### Portfolio Greeks

```q
portfolioGreeks:.portfolio.calculatePortfolioGreeks[tradeTable;marketData;model;config]
```

Returns delta/gamma/theta/vega/rho for European vanilla trades. American and barrier trades return status:\`UNSUPPORTED with null Greeks.

### Portfolio Scenario Risk

```q
portfolioScenarios:.portfolio.generatePortfolioScenarioReport[tradeTable;marketData;model;config]
```

Returns 9 scenarios per trade (36 rows for 4 trades): base, spot +/-1%/+/-5%, vol +/-1pp, rate +/-25bp. Each row includes unitPnL and notionalPnL.

### Scenario Summary

```q
summary:.portfolio.summarizePortfolioRisk[portfolioScenarios]
```

Aggregates notional price and PnL by scenario across all trades.

## Public API

| Function | Returns |
|----------|---------|
| `.engine.priceOption[trade;mkt;model;cfg]` | Price dictionary |
| `.engine.priceOptionWithGrid[trade;mkt;model;cfg]` | Price + full grid |
| `.greeks.calculateGreeks[trade;mkt;model;cfg]` | Greeks table |
| `.risk.generateScenarioReport[trade;mkt;model;cfg]` | Scenario risk |
| `.portfolio.priceTrades[tbl;mkt;model;cfg]` | Portfolio prices |
| `.portfolio.calculatePortfolioGreeks[tbl;mkt;model;cfg]` | Portfolio Greeks |
| `.portfolio.generatePortfolioScenarioReport[tbl;mkt;model;cfg]` | Portfolio scenarios |
| `.portfolio.summarizePortfolioRisk[scenarioTbl]` | Scenario summary |
| `.american.extractEarlyExerciseBoundary[trade;mkt;model;cfg]` | Exercise boundary |
| `.validation.validateEuropeanOption[trade;mkt;model;cfg]` | FDM vs BS price |

## Architecture

```
lib/
  init.q          Silent loader
  utilities.q     Validation, interpolation
  config.q        FDM configuration (explicit + CN)
  product.q       Trade definition (vanilla + barrier)
  market.q        Market data (flat + local vol)
  model.q         BS + local vol model
  grid.q          Grid construction
  payoff.q        Terminal payoff
  boundary.q      European boundaries + barrier
  solver.q        Explicit (flat+localvol) + CN + Thomas
  engine.q        Pricing engine (method routing)
  greeks.q        Greeks
  validation.q    BS closed form + validation
  risk.q          Scenario risk
  american.q      American put analysis
  portfolio.q     Portfolio batch pricing and risk
```

## Running

```
q examples/smoke_test_european_call.q         # European call
q examples/price_portfolio.q                  # Portfolio pricing
q examples/compare_explicit_crank_nicolson.q  # Explicit vs CN
q examples/price_american_put.q               # American put
q examples/price_barrier_options.q            # Barrier options
q examples/price_local_volatility.q           # Local vol
q examples/price_local_volatility_skew.q      # Local vol skew
q tests/run_all_tests.q                       # full suite (24 tests)
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
| v0.7 | Local volatility for European vanilla (explicit) |
| v0.7.1 | Local vol skew sanity + time dependence tests |
| v0.8 | Portfolio pricing and batch risk |

## Roadmap

| Version | Planned |
|---------|---------|
| v0.9 | Multi-underlying portfolios or real market data |
| v1.0 | Production hardening, HDB integration |
