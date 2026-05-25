# qFDM: kdb+/q Finite Difference Pricing Framework

qFDM is a modular kdb+/q pricing framework for equity options using finite-difference methods.

## Current Scope (v0.7.1)

**Supported:**
- European call and put (explicit + Crank-Nicolson)
- American put (explicit only)
- Up-and-out call, down-and-out put (explicit only)
- Local volatility sigma(S,t) for European vanilla (explicit only)
- Greeks, scenario risk, exercise boundary extraction
- Black-Scholes analytical validation

**Not yet supported:**
- CN/American/barrier with local vol
- Dupire calibration, implied vol surface
- Implicit method, portfolio pricing

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

## Local Volatility (v0.7)

The solver supports sigma(S,t) at each grid point via a user-supplied q function.

### Usage

```q
localVolFn:{[spotValue;timePoint] 0.2};
localVolMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;localVolFn];
lvModel:.model.createLocalVolatilityModel[];
show .engine.priceOption[trade;localVolMkt;lvModel;config]
```

### Flat Equivalence

If the local vol function returns a constant, the price matches the flat Black-Scholes explicit FDM price exactly.

### Sanity Checks (v0.7.1)

v0.7.1 adds non-flat local volatility behaviour tests proving the solver correctly uses sigma(S,t):

**Downside-only skew** adds vol below strike while keeping flat vol at and above strike:

```q
downsideOnlySkewFn:{[spotValue;timePoint]
    extraVol:0f | 0.001 * (100f - spotValue);
    0.2 + extraVol
 };
```

Since vol is >= 0.2 everywhere (and strictly higher for S < 100), both call and put prices increase. The put increases more because the extra vol is concentrated in the put exercise region.

**Time dependence** confirms the solver passes the time coordinate:

```q
timeDependentFn:{[spotValue;timePoint] 0.15 + 0.10 * timePoint};
```

### Limitations

- Explicit FDM only (no CN)
- European vanilla only
- No Dupire calibration or implied vol surface
- Function must be vectorized over spot

## Public API

| Function | Returns |
|----------|---------|
| `.engine.priceOption[trade;mkt;model;cfg]` | Price dictionary |
| `.engine.priceOptionWithGrid[trade;mkt;model;cfg]` | Price + full grid |
| `.greeks.calculateGreeks[trade;mkt;model;cfg]` | Greeks table |
| `.risk.generateScenarioReport[trade;mkt;model;cfg]` | Scenario risk |
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
```

## Running

```
q examples/smoke_test_european_call.q         # European call
q examples/compare_explicit_crank_nicolson.q  # Explicit vs CN
q examples/price_american_put.q               # American put
q examples/price_barrier_options.q            # Barrier options
q examples/price_local_volatility.q           # Local vol flat equivalence
q examples/price_local_volatility_skew.q      # Local vol with skew
q tests/run_all_tests.q                       # full suite (21 tests)
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

## Roadmap

| Version | Planned |
|---------|---------|
| v0.8 | Benchmarking and memory estimation |
| v0.9 | Portfolio-level pricing or real market data |
