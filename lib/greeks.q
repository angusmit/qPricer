/ greeks.q — option Greeks from FDM grid and bump-and-reprice
/ Conventions:
/   delta: per 1 unit spot move
/   gamma: per 1 unit spot move squared
/   theta: annual (dV/dt, usually negative for long vanilla)
/   vega:  per 1.00 absolute volatility (vega*0.01 = price change for +1pp vol)
/   rho:   per 1.00 absolute rate (rho*0.0001 = price change for +1bp rate)

.greeks.__volBump:0.005;
.greeks.__rateBump:0.0001;

/ --- Public ---

.greeks.calculateGreeks:{[trade;marketData;model;config]
    gridResult:.engine.priceOptionWithGrid[trade;marketData;model;config];
    solverResult:gridResult`solverResult;
    spotGrid:solverResult`spotGrid;
    timeGrid:solverResult`timeGrid;
    valueGrid:solverResult`valueGrid;
    valueAtTimeZero:solverResult`valueAtTimeZero;
    spot:.market.getSpot[marketData;trade`underlying];
    delta:.greeks.__calculateDeltaFromGrid[spot;spotGrid;valueAtTimeZero];
    gamma:.greeks.__calculateGammaFromGrid[spot;spotGrid;valueAtTimeZero];
    theta:.greeks.__calculateThetaFromGrid[spot;spotGrid;timeGrid;valueGrid];
    vega:.greeks.__calculateVegaByBumpAndReprice[trade;marketData;model;config;.greeks.__volBump];
    rho:.greeks.__calculateRhoByBumpAndReprice[trade;marketData;model;config;.greeks.__rateBump];
    ([] tradeId:enlist trade`tradeId; underlying:enlist trade`underlying;
        optionType:enlist trade`optionType;
        delta:enlist delta; gamma:enlist gamma; theta:enlist theta;
        vega:enlist vega; rho:enlist rho)
 };

/ --- Grid Greeks ---

/ Delta via central difference: (V[i+1] - V[i-1]) / (S[i+1] - S[i-1])
.greeks.__calculateDeltaFromGrid:{[spot;spotGrid;valueAtTimeZero]
    spotIndex:.greeks.__findInteriorIndex[spot;spotGrid];
    (valueAtTimeZero[spotIndex+1] - valueAtTimeZero[spotIndex-1]) % spotGrid[spotIndex+1] - spotGrid[spotIndex-1]
 };

/ Gamma via central second difference: (V[i+1] - 2V[i] + V[i-1]) / dS^2
.greeks.__calculateGammaFromGrid:{[spot;spotGrid;valueAtTimeZero]
    spotIndex:.greeks.__findInteriorIndex[spot;spotGrid];
    dS:spotGrid[spotIndex+1] - spotGrid spotIndex;
    (valueAtTimeZero[spotIndex+1] + valueAtTimeZero[spotIndex-1] - 2f * valueAtTimeZero spotIndex) % dS * dS
 };

/ Theta: (V[spot,t1] - V[spot,t0]) / (t1 - t0), annual
.greeks.__calculateThetaFromGrid:{[spot;spotGrid;timeGrid;valueGrid]
    dt:timeGrid[1] - timeGrid 0;
    priceT0:.utilities.linearInterpolate[spot;spotGrid;valueGrid[;0]];
    priceT1:.utilities.linearInterpolate[spot;spotGrid;valueGrid[;1]];
    (priceT1 - priceT0) % dt
 };

/ Find nearest interior grid index; error if at boundary
.greeks.__findInteriorIndex:{[spot;spotGrid]
    idx:.utilities.findNearestIndex[spot;spotGrid];
    idx:1 | idx;
    idx:idx & (-2) + count spotGrid;
    idx
 };

/ --- Bump-and-reprice Greeks ---

/ Vega via central vol bump; returns per 1.00 absolute vol
.greeks.__calculateVegaByBumpAndReprice:{[trade;marketData;model;config;halfBump]
    mktUp:.market.bumpVolatility[marketData;halfBump];
    mktDn:.market.bumpVolatility[marketData;neg halfBump];
    pUp:.greeks.__getUnitPrice[.engine.priceOption[trade;mktUp;model;config]];
    pDn:.greeks.__getUnitPrice[.engine.priceOption[trade;mktDn;model;config]];
    (pUp - pDn) % 2f * halfBump
 };

/ Rho via central rate bump; returns per 1.00 absolute rate
.greeks.__calculateRhoByBumpAndReprice:{[trade;marketData;model;config;halfBump]
    mktUp:.market.bumpRiskFreeRate[marketData;halfBump];
    mktDn:.market.bumpRiskFreeRate[marketData;neg halfBump];
    pUp:.greeks.__getUnitPrice[.engine.priceOption[trade;mktUp;model;config]];
    pDn:.greeks.__getUnitPrice[.engine.priceOption[trade;mktDn;model;config]];
    (pUp - pDn) % 2f * halfBump
 };

.greeks.__getUnitPrice:{[priceResult] priceResult`unitPrice};
