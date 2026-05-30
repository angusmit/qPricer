/ heston.q - Heston stochastic volatility Monte Carlo (v0.20)

/ --- Validation ---

.heston.validateHestonParams:{[hestonParams]
    if[hestonParams[`initialVariance]<0f; '"initialVariance must be non-negative"];
    if[not hestonParams[`longRunVariance]>0f; '"longRunVariance must be positive"];
    if[not hestonParams[`meanReversion]>0f; '"meanReversion must be positive"];
    if[hestonParams[`volOfVol]<0f; '"volOfVol must be non-negative"];
    corrVal:hestonParams`correlation;
    if[(corrVal< -1f) or corrVal>1f; '"correlation must be between -1 and 1"];
 };

/ --- Correlated Brownian motions ---

.heston.generateCorrelatedBrownian:{[pathCountValue;timeStepCountValue;rho;randomSeedValue]
    totalNormals:pathCountValue*timeStepCountValue;
    allNormals:.montecarlo.__generateNormals[2*totalNormals;randomSeedValue];
    z1Flat:totalNormals#allNormals;
    z2IndFlat:totalNormals _ allNormals;
    z1Matrix:timeStepCountValue cut z1Flat;
    z2IndMatrix:timeStepCountValue cut z2IndFlat;
    rhoComplement:sqrt 1f-rho*rho;
    z2Matrix:(rho*z1Matrix)+rhoComplement*z2IndMatrix;
    `z1`z2!(z1Matrix;z2Matrix)
 };

/ --- Heston path simulation (full truncation Euler) ---

.heston.simulateHestonPaths:{[spotVal;hestonParams;expiry;mcConfig]
    .heston.validateHestonParams hestonParams;
    pathCountValue:mcConfig`pathCount;
    timeStepCountValue:mcConfig`timeStepCount;
    randomSeedValue:mcConfig`randomSeed;
    kappa:hestonParams`meanReversion;
    theta:hestonParams`longRunVariance;
    xi:hestonParams`volOfVol;
    rho:hestonParams`correlation;
    riskFreeRate:hestonParams`riskFreeRate;
    dividendYield:hestonParams`dividendYield;
    v0:hestonParams`initialVariance;
    dtVal:expiry%timeStepCountValue;
    sqrtDt:sqrt dtVal;
    brownianPair:.heston.generateCorrelatedBrownian[pathCountValue;timeStepCountValue;rho;randomSeedValue];
    z1Matrix:brownianPair`z1;
    z2Matrix:brownianPair`z2;
    currentSpot:pathCountValue#spotVal;
    currentVar:pathCountValue#v0;
    spotPathList:enlist currentSpot;
    varPathList:enlist currentVar;
    negVarTotal:0;
    stepIdx:0;
    while[stepIdx<timeStepCountValue;
        z1Step:z1Matrix[;stepIdx];
        z2Step:z2Matrix[;stepIdx];
        vPositive:0f|currentVar;
        sqrtV:sqrt vPositive;
        halfVPos:0.5*vPositive;
        spotDrift:((riskFreeRate-dividendYield)-halfVPos)*dtVal;
        spotDiffusion:sqrtV*sqrtDt*z1Step;
        currentSpot:currentSpot*exp spotDrift+spotDiffusion;
        varDrift:kappa*(theta-vPositive)*dtVal;
        varDiffusion:xi*sqrtV*sqrtDt*z2Step;
        currentVar:currentVar+varDrift+varDiffusion;
        negVarTotal+:sum currentVar<0f;
        spotPathList:spotPathList,enlist currentSpot;
        varPathList:varPathList,enlist currentVar;
        stepIdx+:1];
    spotPaths:flip spotPathList;
    variancePaths:flip varPathList;
    `spotPaths`variancePaths`finalSpot`finalVariance`negativeRawVarianceCount`status`errorMessage!(
        spotPaths;variancePaths;currentSpot;currentVar;negVarTotal;`OK;"")
 };

/ --- European payoff ---

.heston.europeanPayoff:{[spotFinalVector;strike;optionType]
    if[optionType~`call; :0f|spotFinalVector-strike];
    if[optionType~`put; :0f|strike-spotFinalVector];
    '"Unsupported optionType: ",string optionType
 };

/ --- European pricing ---

.heston.priceEuropean:{[trade;marketData;configDict]
    hestonParams:configDict`hestonParams;
    mcConfig:$[`mcConfig in key configDict; configDict`mcConfig; .montecarlo.defaultMcConfig[]];
    .montecarlo.validateMcConfig mcConfig;
    spotVal:marketData`spot;
    expiry:trade`expiry;
    pathResult:.heston.simulateHestonPaths[spotVal;hestonParams;expiry;mcConfig];
    if[not pathResult[`status]~`OK; :pathResult];
    payoffVector:.heston.europeanPayoff[pathResult`finalSpot;trade`strike;trade`optionType];
    riskFreeRate:hestonParams`riskFreeRate;
    priceResult:.montecarlo.priceFromPayoffs[payoffVector;riskFreeRate;expiry;mcConfig`confidenceLevel];
    `tradeId`underlying`optionType`productType`modelType`unitPrice`notionalPrice`method`modelName`standardError`lowerConfidence`upperConfidence`pathCount`status`statusMessage!(
        trade`tradeId;trade`underlying;trade`optionType;`equityOption;`heston;
        priceResult`price;priceResult[`price]*trade`notional;
        `monteCarlo;`heston;
        priceResult`standardError;priceResult`lowerConfidence;priceResult`upperConfidence;
        mcConfig`pathCount;`OK;"")
 };

/ --- Greeks ---

.heston.bumpGreek:{[trade;marketData;configDict;greekName]
    baseResult:.heston.priceEuropean[trade;marketData;configDict];
    basePrice:baseResult`unitPrice;
    bumpedConfig:configDict;
    bumpedMkt:marketData;
    bumpSize:0f;
    hestonParams:configDict`hestonParams;
    if[greekName~`delta;
        bumpSize:marketData[`spot]*0.01;
        bumpedMkt:@[marketData;`spot;:;marketData[`spot]+bumpSize]];
    if[greekName~`vega;
        bumpSize:0.01;
        currentSigma:sqrt hestonParams`initialVariance;
        bumpedSigma:currentSigma+bumpSize;
        bumpedParams:@[hestonParams;`initialVariance;:;bumpedSigma*bumpedSigma];
        bumpedConfig:@[configDict;`hestonParams;:;bumpedParams]];
    if[greekName~`rho;
        bumpSize:0.0001;
        bumpedParams:@[hestonParams;`riskFreeRate;:;hestonParams[`riskFreeRate]+bumpSize];
        bumpedConfig:@[configDict;`hestonParams;:;bumpedParams]];
    if[bumpSize=0f; '"Unsupported greekName: ",string greekName];
    bumpedResult:.heston.priceEuropean[trade;bumpedMkt;bumpedConfig];
    (bumpedResult[`unitPrice]-basePrice)%bumpSize
 };

/ --- Path diagnostics ---

.heston.pathDiagnostics:{[pathResult]
    finalSpotVector:pathResult`finalSpot;
    finalVarVector:pathResult`finalVariance;
    `pathCount`averageFinalSpot`averageFinalVariance`minFinalVariance`maxFinalVariance`negativeRawVarianceCount!(
        count finalSpotVector;
        avg finalSpotVector;
        avg finalVarVector;
        min finalVarVector;
        max finalVarVector;
        pathResult`negativeRawVarianceCount)
 };
