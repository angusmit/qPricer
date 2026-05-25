/ solver.q - FDM solvers for the Black-Scholes PDE
/ Explicit: forward-in-time, supports European/American/barrier
/ Crank-Nicolson: unconditionally stable, European vanilla only (v0.6)

/ ============================================================================
/ Explicit solver (unchanged from v0.5)
/ ============================================================================

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
    currentValues:.payoff.calculateIntrinsicValue[trade;sGrid];
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
        if[hasBarrier; currentValues:.boundary.applyBarrierCondition[trade;sGrid;currentValues]];
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

/ ============================================================================
/ Crank-Nicolson solver (v0.6 - European vanilla only)
/ ============================================================================

.solver.solveCrankNicolson:{[trade;marketData;model;config]
    .solver.__validateCrankNicolsonInputs[trade;marketData;model;config];
    solverInputs:.solver.__buildSolverInputs[trade;marketData;model;config];
    gridDict:solverInputs`grid;
    nTime:config`numberOfTimeSteps;
    expiry:trade`expiry;
    tGrid:gridDict`timeGrid;
    sGrid:gridDict`spotGrid;
    spotStep:gridDict`spotStep;
    timeStep:gridDict`timeStep;
    / Build CN coefficients (constant for flat BS)
    cnCoeffs:.solver.__buildCrankNicolsonCoefficients[solverInputs];
    / Terminal payoff
    currentValues:.payoff.calculateIntrinsicValue[trade;sGrid];
    doFullGrid:0b;
    if[config`returnFullGrid; doFullGrid:1b];
    accumulator:enlist currentValues;
    / Step backward from expiry to t=0
    stepIdx:nTime-1;
    while[stepIdx>=0;
        remainingTime:expiry-tGrid stepIdx;
        currentValues:.solver.__stepCrankNicolson[trade;marketData;gridDict;currentValues;cnCoeffs;remainingTime];
        if[doFullGrid; accumulator:accumulator,enlist currentValues];
        stepIdx-:1];
    valueGrid:(count sGrid)#enlist currentValues;
    if[doFullGrid; valueGrid:flip reverse accumulator];
    outputMeta:`method`modelName`numberOfSpotSteps`numberOfTimeSteps`spotStep`timeStep`stabilityCheckPassed`exerciseStyle!(
        `crankNicolson; model`modelName; config`numberOfSpotSteps; config`numberOfTimeSteps;
        spotStep; timeStep; 1b; trade`exerciseStyle);
    `spotGrid`timeGrid`valueGrid`valueAtTimeZero`metadata!(sGrid;tGrid;valueGrid;currentValues;outputMeta)
 };

/ ============================================================================
/ Shared helpers
/ ============================================================================

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

/ ============================================================================
/ Explicit step
/ ============================================================================

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

/ ============================================================================
/ Crank-Nicolson helpers
/ ============================================================================

.solver.__validateCrankNicolsonInputs:{[trade;marketData;model;config]
    .product.validateOptionTrade trade;
    .market.validateFlatMarketData marketData;
    .model.validateModel model;
    .config.validateFiniteDifferenceConfig config;
    if[not trade[`exerciseStyle]~`european;
        '"Crank-Nicolson v0.6 only supports European vanilla options"];
    if[.product.isBarrierOption trade;
        '"Crank-Nicolson v0.6 only supports European vanilla options"];
 };

/ Build CN operator coefficients for interior nodes (constant for flat BS)
.solver.__buildCrankNicolsonCoefficients:{[solverInputs]
    vol:solverInputs`volatility;
    rate:solverInputs`riskFreeRate;
    divY:solverInputs`dividendYield;
    dS:solverInputs`spotStep;
    dt:solverInputs`timeStep;
    interiorSpot:solverInputs`interiorSpotGrid;
    volSq:vol*vol;
    spotSq:interiorSpot*interiorSpot;
    dsSq:dS*dS;
    driftRate:rate-divY;
    / Spatial operator coefficients per interior node
    diffusionCoeff:0.5*volSq*spotSq%dsSq;
    driftCoeff:driftRate*interiorSpot%(2f*dS);
    lowerOp:diffusionCoeff-driftCoeff;
    mainOp:neg[2f*diffusionCoeff]-rate;
    upperOp:diffusionCoeff+driftCoeff;
    halfDt:0.5*dt;
    / A matrix coefficients (implicit side): I - 0.5*dt*L
    aLowerFull:neg halfDt*lowerOp;
    aMainFull:1f-halfDt*mainOp;
    aUpperFull:neg halfDt*upperOp;
    / B coefficients (explicit side): I + 0.5*dt*L
    bLowerFull:halfDt*lowerOp;
    bMainFull:1f+halfDt*mainOp;
    bUpperFull:halfDt*upperOp;
    / Extract tridiagonal diagonals for A
    triLower:1_aLowerFull;
    triUpper:(-1)_aUpperFull;
    `aLowerFull`aMainFull`aUpperFull`bLowerFull`bMainFull`bUpperFull`triLower`triMain`triUpper!(
        aLowerFull;aMainFull;aUpperFull;bLowerFull;bMainFull;bUpperFull;triLower;aMainFull;triUpper)
 };

/ One CN backward step: solve for previous-time values given next-time values
.solver.__stepCrankNicolson:{[trade;marketData;gridDict;nextFullValues;cnCoeffs;remainingTime]
    sGrid:gridDict`spotGrid;
    interiorNodeCount:(-2)+count sGrid;
    / Compute boundary values at the previous (earlier) time
    lowBoundary:.solver.__cnLowBoundary[trade;marketData;remainingTime];
    highBoundary:.solver.__cnHighBoundary[trade;marketData;gridDict;remainingTime];
    / Build RHS from next (later) time values
    nextBelow:(-2)_nextFullValues;
    nextAt:(-1)_1_nextFullValues;
    nextAbove:2_nextFullValues;
    rhsVector:(cnCoeffs[`bLowerFull]*nextBelow)+(cnCoeffs[`bMainFull]*nextAt)+cnCoeffs[`bUpperFull]*nextAbove;
    / Adjust RHS for known boundaries at previous time
    rhsVector:@[rhsVector;0;-;cnCoeffs[`aLowerFull][0]*lowBoundary];
    rhsVector:@[rhsVector;interiorNodeCount-1;-;cnCoeffs[`aUpperFull][interiorNodeCount-1]*highBoundary];
    / Solve tridiagonal system
    interiorSolution:.solver.__solveTridiagonalSystem[cnCoeffs`triLower;cnCoeffs`triMain;cnCoeffs`triUpper;rhsVector];
    / Assemble full vector
    (enlist lowBoundary),interiorSolution,enlist highBoundary
 };

/ European boundary values for CN
.solver.__cnLowBoundary:{[trade;marketData;remainingTime]
    rDisc:exp neg marketData[`riskFreeRate]*remainingTime;
    if[trade[`optionType]~`call; :0f];
    trade[`strike]*rDisc
 };

.solver.__cnHighBoundary:{[trade;marketData;gridDict;remainingTime]
    rDisc:exp neg marketData[`riskFreeRate]*remainingTime;
    qDisc:exp neg marketData[`dividendYield]*remainingTime;
    sMax:last gridDict`spotGrid;
    if[trade[`optionType]~`call; :(sMax*qDisc)-trade[`strike]*rDisc];
    0f
 };

/ ============================================================================
/ Thomas algorithm tridiagonal solver
/ ============================================================================
/ Solves Ax = d where A is tridiagonal.
/ lowerDiag: length m-1 (sub-diagonal)
/ mainDiag:  length m   (main diagonal)
/ upperDiag: length m-1 (super-diagonal)
/ rhsVector: length m   (right-hand side)
/ Returns:   length m   (solution)

.solver.__solveTridiagonalSystem:{[lowerDiag;mainDiag;upperDiag;rhsVector]
    systemSize:count mainDiag;
    / Forward sweep
    modUpper:(systemSize-1)#0f;
    modRhs:systemSize#0f;
    modUpper[0]:upperDiag[0]%mainDiag[0];
    modRhs[0]:rhsVector[0]%mainDiag[0];
    sweepIdx:1;
    while[sweepIdx<systemSize;
        denom:mainDiag[sweepIdx]-lowerDiag[sweepIdx-1]*modUpper[sweepIdx-1];
        if[sweepIdx<systemSize-1; modUpper[sweepIdx]:upperDiag[sweepIdx]%denom];
        modRhs[sweepIdx]:(rhsVector[sweepIdx]-lowerDiag[sweepIdx-1]*modRhs[sweepIdx-1])%denom;
        sweepIdx+:1];
    / Back substitution
    solution:systemSize#0f;
    solution[systemSize-1]:modRhs[systemSize-1];
    backIdx:systemSize-2;
    while[backIdx>=0;
        solution[backIdx]:modRhs[backIdx]-modUpper[backIdx]*solution[backIdx+1];
        backIdx-:1];
    solution
 };
