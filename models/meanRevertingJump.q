/ meanRevertingJump.q - mean-reverting log-price model with compound Poisson jumps (v0.35)
/ Risk-neutral dynamics:
/     d X_t = kappa * (theta - X_t) dt + sigma dW_t + J_t dN_t
/     S_t = exp X_t
/ N_t is Poisson with intensity lambda; J_t ~ N(jumpMean, jumpVolatility^2).
/ Useful as a foundation for power/electricity spike modelling and oil-shock scenarios.
/ v0.35 is simulation-only: no closed-form pricing, no calibration, no real data.
/ All parameters are interpreted as risk-neutral; market price of risk is deferred.
/ Initial state x0 is a function argument, matching v0.33/v0.34 conventions.
/ Public:  validateParams, transitionMomentsApprox, simulatePaths, pathDiagnostics,
/          jumpDiagnostics, europeanOptionPriceMC, spikeProbability, scenarioPaths,
/          meanPath, expectedJumpCount.
/ Private: __exactBaseStep, __effectiveBaseVariance, __jumpCounts, __jumpSizes,
/          __stepWithJumps, __payoff.

/ Threshold below which |kappa * tau| triggers Taylor branches.
.commodity.mrjump.__smallKappaTau:1e-6;

.commodity.mrjump.validateParams:{[params]
    requiredKeys:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility;
    paramKeys:key params;
    missingKeys:requiredKeys where not requiredKeys in paramKeys;
    if[0<count missingKeys;
        missingStr:", " sv string missingKeys;
        '"mrjump params missing: ",missingStr];
    if[params[`meanReversionSpeed]<=0f; '"mrjump meanReversionSpeed must be positive"];
    if[params[`volatility]<0f; '"mrjump volatility must be non-negative"];
    if[params[`jumpIntensity]<0f; '"mrjump jumpIntensity must be non-negative"];
    if[params[`jumpVolatility]<0f; '"mrjump jumpVolatility must be non-negative"];
 };

/ Var of OU base over tau: sigma^2 * (1 - exp(-2 kappa tau)) / (2 kappa).
/ Small-kappa Taylor: sigma^2 * tau * (1 - kt + (2/3) kt^2).
.commodity.mrjump.__effectiveBaseVariance:{[params;tau]
    kappaVal:params`meanReversionSpeed;
    sigmaVal:params`volatility;
    sigmaSq:sigmaVal*sigmaVal;
    kappaTau:kappaVal*tau;
    if[(abs kappaTau)<.commodity.mrjump.__smallKappaTau;
        taylorFactor:(1f-kappaTau)+(2f%3f)*kappaTau*kappaTau;
        :sigmaSq*tau*taylorFactor];
    sigmaSq*(1f-exp neg 2f*kappaTau)%(2f*kappaVal)
 };

/ Conditional log-price moments under compound Poisson + OU.
/ m = exp(-kt) X0 + (1-exp(-kt)) theta + lambda * jumpMean * (1-exp(-kt))/k
/ v = (sigma^2 + lambda * (jumpMean^2 + jumpVol^2)) * (1-exp(-2kt))/(2k)
/ Taylor branches for |kt| < smallKappaTau.
.commodity.mrjump.transitionMomentsApprox:{[x0;params;tau]
    .commodity.mrjump.validateParams params;
    if[tau<0f; '"mrjump tau must be non-negative"];
    kappaVal:params`meanReversionSpeed;
    thetaVal:params`longRunLogMean;
    sigmaVal:params`volatility;
    lambdaVal:params`jumpIntensity;
    jumpMeanVal:params`jumpMean;
    jumpVolVal:params`jumpVolatility;
    kappaTau:kappaVal*tau;
    smallKt:(abs kappaTau)<.commodity.mrjump.__smallKappaTau;
    decayVal:exp neg kappaVal*tau;
    decay2Val:exp neg 2f*kappaVal*tau;
    factor1:$[smallKt;
        tau*(1f-0.5*kappaTau)+(kappaTau*kappaTau)%6f;
        (1f-decayVal)%kappaVal];
    factor2:$[smallKt;
        tau*(1f-kappaTau)+(2f%3f)*kappaTau*kappaTau;
        (1f-decay2Val)%(2f*kappaVal)];
    meanVal:(decayVal*x0)+((1f-decayVal)*thetaVal)+lambdaVal*jumpMeanVal*factor1;
    varianceVal:((sigmaVal*sigmaVal)+lambdaVal*(jumpMeanVal*jumpMeanVal)+jumpVolVal*jumpVolVal)*factor2;
    `mean`variance!(meanVal;varianceVal)
 };

/ Expected jump count over [0, expiry] under Poisson process.
.commodity.mrjump.expectedJumpCount:{[params;expiry]
    (params`jumpIntensity)*expiry
 };

/ Probability of at least one jump before expiry.
.commodity.mrjump.spikeProbability:{[params;expiry]
    1f-exp neg (params`jumpIntensity)*expiry
 };

/ Mean log-price path E[X_tau] across a tenor vector (returns log prices, not spot).
.commodity.mrjump.meanPath:{[x0;params;tenorVector]
    .commodity.mrjump.validateParams params;
    {[xInit;p;t] (.commodity.mrjump.transitionMomentsApprox[xInit;p;t])`mean}[x0;params;] each tenorVector
 };

/ Exact OU step for the diffusion component (no jumps yet).
.commodity.mrjump.__exactBaseStep:{[currentLogPrice;params;dtVal;normalDraws]
    kappaVal:params`meanReversionSpeed;
    thetaVal:params`longRunLogMean;
    decayVal:exp neg kappaVal*dtVal;
    stepVarianceVal:.commodity.mrjump.__effectiveBaseVariance[params;dtVal];
    stepStdVal:sqrt stepVarianceVal;
    (decayVal*currentLogPrice)+((1f-decayVal)*thetaVal)+stepStdVal*normalDraws
 };

/ Poisson jump count matrix (pathCount x stepCount) via inverse-CDF method.
/ Reuses .merton.__poissonFromUniforms which is capped at 19 terms (safe up to lambda*dt ~ 10).
.commodity.mrjump.__jumpCounts:{[pathCount;stepCount;jumpIntensity;dtVal;randomSeed]
    lambdaDt:jumpIntensity*dtVal;
    totalCells:pathCount*stepCount;
    system "S ",string randomSeed;
    uniforms:totalCells?1f;
    counts:.merton.__poissonFromUniforms[lambdaDt;uniforms];
    stepCount cut counts
 };

/ Aggregate jump-size matrix. For K_ij ~ Poisson, sum of K_ij draws of N(jumpMean, jumpVol^2)
/ is exactly N(K*jumpMean, K*jumpVol^2), so we sample with K*jumpMean + sqrt(K)*jumpVol*Z.
.commodity.mrjump.__jumpSizes:{[jumpCountMatrix;params;randomSeed]
    jumpMeanVal:params`jumpMean;
    jumpVolVal:params`jumpVolatility;
    pathCountValue:count jumpCountMatrix;
    stepCountValue:count jumpCountMatrix 0;
    totalCells:pathCountValue*stepCountValue;
    flatNormals:.montecarlo.__generateNormals[totalCells;randomSeed];
    normalMatrix:stepCountValue cut flatNormals;
    (jumpMeanVal*jumpCountMatrix)+jumpVolVal*(sqrt jumpCountMatrix)*normalMatrix
 };

/ One full step: OU base plus aggregated jump for the timestep.
.commodity.mrjump.__stepWithJumps:{[currentLogPrice;params;dtVal;normalDraws;jumpCounts;jumpSizes]
    baseAfter:.commodity.mrjump.__exactBaseStep[currentLogPrice;params;dtVal;normalDraws];
    baseAfter+jumpSizes
 };

.commodity.mrjump.__payoff:{[optType;terminalPrice;strikePrice]
    $[optType=`call;0f|terminalPrice-strikePrice;0f|strikePrice-terminalPrice]
 };

/ Simulate compound-Poisson + OU paths. Returns dict:
/   pathMatrix     - pathCount x stepCount spot prices (excluding t=0)
/   jumpCountMatrix - pathCount x stepCount integer jump counts per step
/   jumpSizeMatrix  - pathCount x stepCount aggregated jump contributions per step
.commodity.mrjump.simulatePaths:{[x0;params;expiry;mcConfig]
    .commodity.mrjump.validateParams params;
    .montecarlo.validateMcConfig mcConfig;
    if[expiry<=0f; '"mrjump expiry must be positive"];
    pathCountValue:mcConfig`pathCount;
    stepCountValue:mcConfig`timeStepCount;
    seedValue:mcConfig`randomSeed;
    useAntithetic:mcConfig`antithetic;
    dtVal:expiry%stepCountValue;
    diffusionNormals:.montecarlo.generateNormalPaths[pathCountValue;stepCountValue;seedValue];
    if[useAntithetic; diffusionNormals:.montecarlo.generateAntitheticNormals[pathCountValue;stepCountValue;seedValue]];
    jumpCountMatrix:.commodity.mrjump.__jumpCounts[pathCountValue;stepCountValue;params`jumpIntensity;dtVal;seedValue+1];
    jumpSizeMatrix:.commodity.mrjump.__jumpSizes[jumpCountMatrix;params;seedValue+2];
    currentLogPrices:pathCountValue#x0;
    logPathRows:enlist currentLogPrices;
    stepIdx:0;
    while[stepIdx<stepCountValue;
        normalsAtStep:diffusionNormals[;stepIdx];
        jumpCountsAtStep:jumpCountMatrix[;stepIdx];
        jumpSizesAtStep:jumpSizeMatrix[;stepIdx];
        currentLogPrices:.commodity.mrjump.__stepWithJumps[currentLogPrices;params;dtVal;normalsAtStep;jumpCountsAtStep;jumpSizesAtStep];
        logPathRows:logPathRows,enlist currentLogPrices;
        stepIdx+:1];
    pathMatrix:flip exp 1_logPathRows;
    `pathMatrix`jumpCountMatrix`jumpSizeMatrix!(pathMatrix;jumpCountMatrix;jumpSizeMatrix)
 };

/ Standard path diagnostics from .pathdiag.summary.
.commodity.mrjump.pathDiagnostics:{[pathMatrix]
    .pathdiag.summary[pathMatrix;1f]
 };

/ Combined price-and-jump diagnostics for a simulatePaths result.
.commodity.mrjump.jumpDiagnostics:{[pathResult]
    pathMatrix:pathResult`pathMatrix;
    jumpCountMatrix:pathResult`jumpCountMatrix;
    pathCountValue:count pathMatrix;
    stepCountValue:count pathMatrix 0;
    finalPrices:last each pathMatrix;
    minPrices:min each pathMatrix;
    maxPrices:max each pathMatrix;
    totalJumpsPerPath:sum each jumpCountMatrix;
    averageJumpCountValue:avg totalJumpsPerPath;
    maxJumpCountValue:max totalJumpsPerPath;
    jumpPathFractionValue:(sum totalJumpsPerPath>0)%pathCountValue;
    `pathCount`timeStepCount`averageFinalPrice`averageMinimumPrice`averageMaximumPrice`averageJumpCount`maxJumpCount`jumpPathFraction`status`errorMessage!(
        pathCountValue;
        stepCountValue;
        avg finalPrices;
        avg minPrices;
        avg maxPrices;
        averageJumpCountValue;
        maxJumpCountValue;
        jumpPathFractionValue;
        `OK;
        "")
 };

/ Monte Carlo European option price.
.commodity.mrjump.europeanOptionPriceMC:{[optType;x0;strikePrice;expiry;riskFreeRate;params;mcConfig]
    if[not optType in `call`put; '"mrjump option type must be call or put"];
    if[strikePrice<=0f; '"mrjump strike must be positive"];
    pathResult:.commodity.mrjump.simulatePaths[x0;params;expiry;mcConfig];
    pathMatrix:pathResult`pathMatrix;
    terminalSpots:last each pathMatrix;
    payoffVector:.commodity.mrjump.__payoff[optType;terminalSpots;strikePrice];
    .montecarlo.priceFromPayoffs[payoffVector;riskFreeRate;expiry;mcConfig`confidenceLevel]
 };

/ Apply simple scenario bumps then resimulate. Recognised scenarioConfig keys:
/   jumpIntensityBump, jumpMeanBump, volatilityBump, longRunLogMeanShift.
/ Unrecognised keys are ignored; returns the same dict shape as simulatePaths.
.commodity.mrjump.scenarioPaths:{[x0;params;expiry;mcConfig;scenarioConfig]
    scenarioParams:params;
    scenarioKeys:key scenarioConfig;
    if[`jumpIntensityBump in scenarioKeys;
        scenarioParams:@[scenarioParams;`jumpIntensity;+;scenarioConfig`jumpIntensityBump]];
    if[`jumpMeanBump in scenarioKeys;
        scenarioParams:@[scenarioParams;`jumpMean;+;scenarioConfig`jumpMeanBump]];
    if[`volatilityBump in scenarioKeys;
        scenarioParams:@[scenarioParams;`volatility;+;scenarioConfig`volatilityBump]];
    if[`longRunLogMeanShift in scenarioKeys;
        scenarioParams:@[scenarioParams;`longRunLogMean;+;scenarioConfig`longRunLogMeanShift]];
    .commodity.mrjump.simulatePaths[x0;scenarioParams;expiry;mcConfig]
 };
