/ solver.q - explicit FDM solver for BS PDE
/ Supports European, American, and barrier options
/ Backward step: V(t) = V(t+dt) + dt*(0.5*vol^2*S^2*gamma + (r-q)*S*delta - r*V)

.solver.solveExplicitFiniteDifference:{[trade;marketData;model;config]
    .product.validateOptionTrade trade;
    .market.validateFlatMarketData marketData;
    .model.validateModel model;
    .config.validateFiniteDifferenceConfig config;
    solverInputs:.solver.__buildSolverInputs[trade;marketData;model;config];
    stable:1b;
    if[config`stabilityCheck; stable:.solver.__checkExplicitStability solverInputs];
    gridDict:solverInputs`grid;
    nTime:config`numberOfTimeSteps;
    expiry:trade`expiry;
    tGrid:gridDict`timeGrid;
    sGrid:gridDict`spotGrid;
    / Terminal payoff
    currentValues:.payoff.calculateIntrinsicValue[trade;sGrid];
    / Apply barrier condition to terminal payoff
    currentValues:.boundary.applyBarrierCondition[trade;sGrid;currentValues];
    isAmerican:trade[`exerciseStyle]~`american;
    intrinsicValues:currentValues;
    hasBarrier:.product.isBarrierOption trade;
    doFullGrid:0b;
    if[config`returnFullGrid; doFullGrid:1b];
    accumulator:enlist currentValues;
    stepIdx:nTime-1;
    while[stepIdx>=0;
        currentValues:.solver.__stepExplicit[currentValues;solverInputs];
        currentValues:.boundary.applyEuropeanBoundary[trade;marketData;gridDict;currentValues;expiry-tGrid stepIdx];
        / Barrier knock-out condition
        if[hasBarrier; currentValues:.boundary.applyBarrierCondition[trade;sGrid;currentValues]];
        / American early exercise
        if[isAmerican; currentValues:currentValues|intrinsicValues];
        if[doFullGrid; accumulator:accumulator,enlist currentValues];
        stepIdx-:1];
    valueGrid:(count sGrid)#enlist currentValues;
    if[doFullGrid; valueGrid:flip reverse accumulator];
    outputMeta:`method`modelName`numberOfSpotSteps`numberOfTimeSteps`spotStep`timeStep`stabilityCheckPassed`exerciseStyle!(
        config`method; model`modelName; config`numberOfSpotSteps; config`numberOfTimeSteps;
        gridDict`spotStep; gridDict`timeStep; stable; trade`exerciseStyle);
    `spotGrid`timeGrid`valueGrid`valueAtTimeZero`metadata!(sGrid;tGrid;valueGrid;currentValues;outputMeta)
 };

.solver.__buildSolverInputs:{[trade;marketData;model;config]
    gridDict:.grid.buildFiniteDifferenceGrid[trade;marketData;config];
    sGrid:gridDict`spotGrid;
    `grid`riskFreeRate`dividendYield`volatility`spotStep`timeStep`interiorSpotGrid!(
        gridDict;
        .market.getRiskFreeRate[marketData;trade`expiry];
        .market.getDividendYield[marketData;trade`underlying;trade`expiry];
        .market.getVolatility[marketData;trade`underlying;trade`strike;trade`expiry];
        gridDict`spotStep; gridDict`timeStep; (-1)_1_sGrid)
 };

.solver.__stepExplicit:{[vals;solverInputs]
    vUp:2_vals;
    vAt:(-1)_1_vals;
    vDn:(-2)_vals;
    dS:solverInputs`spotStep;
    dt:solverInputs`timeStep;
    vol:solverInputs`volatility;
    rate:solverInputs`riskFreeRate;
    divY:solverInputs`dividendYield;
    interiorSpot:solverInputs`interiorSpotGrid;
    delta:(vUp-vDn)%2f*dS;
    gamma:(vUp+vDn-2f*vAt)%(dS*dS);
    diffusion:0.5*vol*vol*interiorSpot*interiorSpot*gamma;
    convection:(rate-divY)*interiorSpot*delta;
    discount:rate*vAt;
    continuation:vAt+dt*(diffusion+convection-discount);
    (enlist 0f),continuation,enlist 0f
 };

.solver.__checkExplicitStability:{[solverInputs]
    vol:solverInputs`volatility;
    dS:solverInputs`spotStep;
    dt:solverInputs`timeStep;
    sMax:last solverInputs[`grid]`spotGrid;
    expiry:last solverInputs[`grid]`timeGrid;
    diffCoeff:solverInputs[`riskFreeRate]+vol*vol*sMax*sMax%(dS*dS);
    lambda:dt*diffCoeff;
    if[lambda>1f;
        minSteps:ceiling expiry*diffCoeff;
        '"Explicit FDM unstable: lambda=",string[lambda],". Increase numberOfTimeSteps to at least ",string minSteps];
    1b
 };
