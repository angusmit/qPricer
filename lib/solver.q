/ solver.q — explicit finite-difference solver for the Black-Scholes PDE
/ Backward step: V(t) = V(t+dt) + dt*(0.5*vol^2*S^2*gamma + (r-q)*S*delta - r*V)

.solver.solveExplicitFiniteDifference:{[trade;marketData;model;config]
    .product.validateOptionTrade trade;
    .market.validateFlatMarketData marketData;
    .model.validateModel model;
    .config.validateFiniteDifferenceConfig config;
    si:.solver.__buildSolverInputs[trade;marketData;model;config];
    stable:1b;
    if[config`stabilityCheck; stable:.solver.__checkExplicitStability si];
    grid:si`grid;
    nTime:config`numberOfTimeSteps;
    expiry:trade`expiry;
    tGrid:grid`timeGrid;
    sGrid:grid`spotGrid;
    cv:.payoff.calculateIntrinsicValue[trade;sGrid];
    doFullGrid:0b;
    if[config`returnFullGrid; doFullGrid:1b];
    acc:enlist cv;
    k:nTime-1;
    while[k>=0;
          cv:.solver.__stepExplicit[cv;si];
          cv:.boundary.applyEuropeanBoundary[trade;marketData;grid;cv;expiry-tGrid k];
          if[doFullGrid; acc:acc,enlist cv];
          k-:1];
    vGrid:(count sGrid)#enlist cv;
    if[doFullGrid; vGrid:flip reverse acc];
    md:`method`modelName`numberOfSpotSteps`numberOfTimeSteps`spotStep`timeStep`stabilityCheckPassed!(
        config`method; model`modelName; config`numberOfSpotSteps; config`numberOfTimeSteps;
        grid`spotStep; grid`timeStep; stable);
    `spotGrid`timeGrid`valueGrid`valueAtTimeZero`metadata!(sGrid;tGrid;vGrid;cv;md)
 };

.solver.__buildSolverInputs:{[trade;marketData;model;config]
    grid:.grid.buildFiniteDifferenceGrid[trade;marketData;config];
    sGrid:grid`spotGrid;
    `grid`riskFreeRate`dividendYield`volatility`spotStep`timeStep`interiorSpotGrid!(
        grid;
        .market.getRiskFreeRate[marketData;trade`expiry];
        .market.getDividendYield[marketData;trade`underlying;trade`expiry];
        .market.getVolatility[marketData;trade`underlying;trade`strike;trade`expiry];
        grid`spotStep; grid`timeStep; (-1)_1_sGrid)
 };

.solver.__stepExplicit:{[vals;si]
    vUp:2_vals;
    vAt:(-1)_1_vals;
    vDn:(-2)_vals;
    dS:si`spotStep;
    dt:si`timeStep;
    vol:si`volatility;
    rate:si`riskFreeRate;
    divY:si`dividendYield;
    iS:si`interiorSpotGrid;
    delta:(vUp-vDn)%2f*dS;
    gamma:(vUp+vDn-2f*vAt)%(dS*dS);
    diffusion:0.5*vol*vol*iS*iS*gamma;
    convection:(rate-divY)*iS*delta;
    discount:rate*vAt;
    cont:vAt+dt*(diffusion+convection-discount);
    (enlist 0f),cont,enlist 0f
 };

.solver.__checkExplicitStability:{[si]
    vol:si`volatility;
    dS:si`spotStep;
    dt:si`timeStep;
    sMax:last si[`grid]`spotGrid;
    lambda:dt*(si[`riskFreeRate]+vol*vol*sMax*sMax%(dS*dS));
    if[lambda>1f;
       '"Explicit FDM unstable: lambda=",string[lambda],". Increase numberOfTimeSteps to at least ",string ceiling 1%(si[`riskFreeRate]+vol*vol*sMax*sMax%(dS*dS))];
    1b
 };
