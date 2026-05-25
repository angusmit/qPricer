/ greeks.q — delta/gamma/theta from grid, vega/rho by bump-and-reprice

.greeks.__volBump:0.01;
.greeks.__rateBump:0.0001;

.greeks.calculateGreeks:{[trade;marketData;model;config]
    cfg:@[config;`returnFullGrid;:;1b];
    solverResult:.solver.solveExplicitFiniteDifference[trade;marketData;model;cfg];
    spot:marketData`spot;
    sGrid:solverResult`spotGrid;
    vt0:solverResult`valueAtTimeZero;
    delta:.greeks.__deltaFromGrid[spot;sGrid;vt0];
    gamma:.greeks.__gammaFromGrid[spot;sGrid;vt0];
    theta:.greeks.__thetaFromGrid[spot;solverResult`timeGrid;sGrid;solverResult`valueGrid];
    vega:.greeks.__vegaBump[trade;marketData;model;cfg;.greeks.__volBump];
    rho:.greeks.__rhoBump[trade;marketData;model;cfg;.greeks.__rateBump];
    ([] tradeId:enlist trade`tradeId; underlying:enlist trade`underlying;
        optionType:enlist trade`optionType;
        delta:enlist delta; gamma:enlist gamma; theta:enlist theta;
        vega:enlist vega; rho:enlist rho)
 };

.greeks.__deltaFromGrid:{[spot;sGrid;vt0]
    dS:sGrid[1]-sGrid 0;
    i:(1|.utilities.findNearestIndex[spot;sGrid])&-2+count sGrid;
    (vt0[i+1]-vt0[i-1])%2f*dS
 };

.greeks.__gammaFromGrid:{[spot;sGrid;vt0]
    dS:sGrid[1]-sGrid 0;
    i:(1|.utilities.findNearestIndex[spot;sGrid])&-2+count sGrid;
    (vt0[i+1]+vt0[i-1]-2f*vt0 i)%(dS*dS)
 };

.greeks.__thetaFromGrid:{[spot;tGrid;sGrid;vGrid]
    dt:tGrid[1]-tGrid 0;
    p0:.utilities.linearInterpolate[spot;sGrid;vGrid[;0]];
    p1:.utilities.linearInterpolate[spot;sGrid;vGrid[;1]];
    (p1-p0)%dt
 };

.greeks.__vegaBump:{[trade;marketData;model;config;bump]
    mktUp:.market.bumpVolatility[marketData;bump];
    srUp:.solver.solveExplicitFiniteDifference[trade;mktUp;model;config];
    pUp:.engine.__interpolatePriceFromGrid[marketData`spot;srUp`spotGrid;srUp`valueAtTimeZero;config`interpolationMethod];
    mktDn:.market.bumpVolatility[marketData;neg bump];
    srDn:.solver.solveExplicitFiniteDifference[trade;mktDn;model;config];
    pDn:.engine.__interpolatePriceFromGrid[marketData`spot;srDn`spotGrid;srDn`valueAtTimeZero;config`interpolationMethod];
    (pUp-pDn)%2f*bump
 };

.greeks.__rhoBump:{[trade;marketData;model;config;bump]
    mktUp:.market.bumpRiskFreeRate[marketData;bump];
    srUp:.solver.solveExplicitFiniteDifference[trade;mktUp;model;config];
    pUp:.engine.__interpolatePriceFromGrid[marketData`spot;srUp`spotGrid;srUp`valueAtTimeZero;config`interpolationMethod];
    mktDn:.market.bumpRiskFreeRate[marketData;neg bump];
    srDn:.solver.solveExplicitFiniteDifference[trade;mktDn;model;config];
    pDn:.engine.__interpolatePriceFromGrid[marketData`spot;srDn`spotGrid;srDn`valueAtTimeZero;config`interpolationMethod];
    (pUp-pDn)%2f*bump
 };
