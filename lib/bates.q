/ bates.q - Bates stochastic volatility with jumps (v0.25)

/ --- Validation ---

.bates.validateParams:{[batesParams]
    if[batesParams[`initialVariance]<0f; '"initialVariance must be non-negative"];
    if[not batesParams[`longRunVariance]>0f; '"longRunVariance must be positive"];
    if[not batesParams[`meanReversion]>0f; '"meanReversion must be positive"];
    if[batesParams[`volOfVol]<0f; '"volOfVol must be non-negative"];
    corrVal:batesParams`correlation;
    if[(corrVal<= -1f) or corrVal>=1f; '"correlation must be between -1 and 1"];
    if[batesParams[`jumpIntensity]<0f; '"jumpIntensity must be non-negative"];
    if[batesParams[`jumpVolatility]<0f; '"jumpVolatility must be non-negative"];
 };

.bates.jumpCompensator:{[batesParams]
    muJ:batesParams`jumpMean;
    sigmaJ:batesParams`jumpVolatility;
    halfSigmaJ2:0.5*sigmaJ*sigmaJ;
    (exp muJ+halfSigmaJ2)-1f
 };

/ --- Path simulation (full truncation Euler + jumps) ---

.bates.simulateBatesPaths:{[spotVal;batesParams;expiry;mcConfig]
    .bates.validateParams batesParams;
    pathCountValue:mcConfig`pathCount;
    timeStepCountValue:mcConfig`timeStepCount;
    randomSeedValue:mcConfig`randomSeed;
    kappa:batesParams`meanReversion;
    theta:batesParams`longRunVariance;
    xi:batesParams`volOfVol;
    rho:batesParams`correlation;
    riskFreeRate:batesParams`riskFreeRate;
    dividendYield:batesParams`dividendYield;
    v0:batesParams`initialVariance;
    lambdaVal:batesParams`jumpIntensity;
    muJVal:batesParams`jumpMean;
    sigmaJVal:batesParams`jumpVolatility;
    dtVal:expiry%timeStepCountValue;
    sqrtDt:sqrt dtVal;
    lambdaDt:lambdaVal*dtVal;
    halfSigmaJ2:0.5*sigmaJVal*sigmaJVal;
    kappaJ:(exp muJVal+halfSigmaJ2)-1f;
    lambdaKappa:lambdaVal*kappaJ;
    rMinusQ:riskFreeRate-dividendYield;
    / Correlated Brownians for spot/vol
    brownianPair:.heston.generateCorrelatedBrownian[pathCountValue;timeStepCountValue;rho;randomSeedValue];
    z1Matrix:brownianPair`z1;
    z2Matrix:brownianPair`z2;
    / Jump size normals
    totalNormals:pathCountValue*timeStepCountValue;
    jumpNormals:.montecarlo.__generateNormals[totalNormals;randomSeedValue+10];
    jumpMatrix:timeStepCountValue cut jumpNormals;
    / Poisson uniforms
    system "S ",string randomSeedValue+20;
    allUniforms:totalNormals?1f;
    uniformMatrix:timeStepCountValue cut allUniforms;
    / Initialize
    currentSpot:pathCountValue#spotVal;
    currentVar:pathCountValue#v0;
    totalJumpCounts:pathCountValue#0;
    negVarTotal:0;
    spotPathList:enlist currentSpot;
    varPathList:enlist currentVar;
    / Time stepping
    stepIdx:0;
    while[stepIdx<timeStepCountValue;
        z1Step:z1Matrix[;stepIdx];
        z2Step:z2Matrix[;stepIdx];
        zJump:jumpMatrix[;stepIdx];
        uPoisson:uniformMatrix[;stepIdx];
        / Heston: full truncation
        vPositive:0f|currentVar;
        sqrtV:sqrt vPositive;
        / Variance update
        varDrift:kappa*(theta-vPositive)*dtVal;
        varDiffusion:xi*sqrtV*sqrtDt*z2Step;
        currentVar:currentVar+varDrift+varDiffusion;
        negVarTotal+:sum currentVar<0f;
        / Jumps
        jumpCounts:.merton.__poissonFromUniforms[lambdaDt;uPoisson];
        totalJumpCounts+:jumpCounts;
        floatJumps:`float$jumpCounts;
        sqrtJumps:sqrt floatJumps;
        jumpSums:(floatJumps*muJVal)+sqrtJumps*sigmaJVal*zJump;
        / Spot update (compensated drift with stochastic vol)
        halfVPos:0.5*vPositive;
        spotDrift:((rMinusQ-lambdaKappa)-halfVPos)*dtVal;
        spotDiffusion:sqrtV*sqrtDt*z1Step;
        logReturn:(spotDrift+spotDiffusion)+jumpSums;
        currentSpot:currentSpot*exp logReturn;
        spotPathList:spotPathList,enlist currentSpot;
        varPathList:varPathList,enlist currentVar;
        stepIdx+:1];
    spotPaths:flip spotPathList;
    variancePaths:flip varPathList;
    `spotPaths`variancePaths`finalSpot`finalVariance`totalJumpCount`negativeRawVarianceCount`status`errorMessage!(
        spotPaths;variancePaths;currentSpot;currentVar;totalJumpCounts;negVarTotal;`OK;"")
 };

/ --- European pricing ---

.bates.europeanPayoff:{[finalSpotVector;strike;optionType]
    if[optionType~`call; :0f|finalSpotVector-strike];
    if[optionType~`put; :0f|strike-finalSpotVector];
    '"Unsupported optionType: ",string optionType
 };

.bates.priceEuropean:{[trade;marketData;configDict]
    batesParams:configDict`batesParams;
    mcConfig:$[`mcConfig in key configDict;configDict`mcConfig;.montecarlo.defaultMcConfig[]];
    spotVal:marketData`spot;
    expiryVal:trade`expiry;
    pathResult:.bates.simulateBatesPaths[spotVal;batesParams;expiryVal;mcConfig];
    if[not pathResult[`status]~`OK; :pathResult];
    payoffVector:.bates.europeanPayoff[pathResult`finalSpot;trade`strike;trade`optionType];
    riskFreeRate:batesParams`riskFreeRate;
    priceResult:.montecarlo.priceFromPayoffs[payoffVector;riskFreeRate;expiryVal;mcConfig`confidenceLevel];
    `tradeId`underlying`optionType`productType`modelType`unitPrice`notionalPrice`method`modelName`standardError`lowerConfidence`upperConfidence`pathCount`status`statusMessage!(
        trade`tradeId;trade`underlying;trade`optionType;`equityOption;`bates;
        priceResult`price;priceResult[`price]*trade`notional;
        `monteCarlo;`bates;
        priceResult`standardError;priceResult`lowerConfidence;priceResult`upperConfidence;
        mcConfig`pathCount;`OK;"")
 };

/ --- Greeks ---

.bates.bumpGreek:{[trade;marketData;configDict;greekName]
    baseResult:.bates.priceEuropean[trade;marketData;configDict];
    basePrice:baseResult`unitPrice;
    bumpedConfig:configDict;
    bumpedMkt:marketData;
    bumpSize:0f;
    batesParams:configDict`batesParams;
    if[greekName~`delta;
        bumpSize:marketData[`spot]*0.01;
        bumpedMkt:@[marketData;`spot;:;marketData[`spot]+bumpSize]];
    if[greekName~`vega;
        bumpSize:0.01;
        currentSigma:sqrt batesParams`initialVariance;
        bumpedSigma:currentSigma+bumpSize;
        bumpedParams:@[batesParams;`initialVariance;:;bumpedSigma*bumpedSigma];
        bumpedConfig:@[configDict;`batesParams;:;bumpedParams]];
    if[greekName~`rho;
        bumpSize:0.0001;
        bumpedParams:@[batesParams;`riskFreeRate;:;batesParams[`riskFreeRate]+bumpSize];
        bumpedConfig:@[configDict;`batesParams;:;bumpedParams]];
    if[greekName~`jumpIntensitySensitivity;
        bumpSize:0.01;
        bumpedParams:@[batesParams;`jumpIntensity;:;batesParams[`jumpIntensity]+bumpSize];
        bumpedConfig:@[configDict;`batesParams;:;bumpedParams]];
    if[greekName~`jumpMeanSensitivity;
        bumpSize:0.01;
        bumpedParams:@[batesParams;`jumpMean;:;batesParams[`jumpMean]+bumpSize];
        bumpedConfig:@[configDict;`batesParams;:;bumpedParams]];
    if[greekName~`jumpVolatilitySensitivity;
        bumpSize:0.01;
        bumpedParams:@[batesParams;`jumpVolatility;:;batesParams[`jumpVolatility]+bumpSize];
        bumpedConfig:@[configDict;`batesParams;:;bumpedParams]];
    if[greekName~`volOfVolSensitivity;
        bumpSize:0.01;
        bumpedParams:@[batesParams;`volOfVol;:;batesParams[`volOfVol]+bumpSize];
        bumpedConfig:@[configDict;`batesParams;:;bumpedParams]];
    if[bumpSize=0f; '"Unsupported greekName: ",string greekName];
    bumpedResult:.bates.priceEuropean[trade;bumpedMkt;bumpedConfig];
    (bumpedResult[`unitPrice]-basePrice)%bumpSize
 };

/ --- Path diagnostics ---

.bates.pathDiagnostics:{[pathResult]
    finalSpotVector:pathResult`finalSpot;
    finalVarVector:pathResult`finalVariance;
    jumpCounts:pathResult`totalJumpCount;
    `pathCount`timeStepCount`averageFinalSpot`averageFinalVariance`averageTotalJumpCount`maxTotalJumpCount`minFinalSpot`maxFinalSpot`negativeRawVarianceCount!(
        count finalSpotVector;
        (count (pathResult`spotPaths) 0)-1;
        avg finalSpotVector;
        avg finalVarVector;
        avg `float$jumpCounts;
        max jumpCounts;
        min finalSpotVector;
        max finalSpotVector;
        pathResult`negativeRawVarianceCount)
 };

/ --- Calibration ---

.bates.parameterGrid:{[paramGridDict]
    ivL:paramGridDict`initialVarianceList;
    lvL:paramGridDict`longRunVarianceList;
    mrL:paramGridDict`meanReversionList;
    vvL:paramGridDict`volOfVolList;
    coL:paramGridDict`correlationList;
    jiL:paramGridDict`jumpIntensityList;
    jmL:paramGridDict`jumpMeanList;
    jvL:paramGridDict`jumpVolatilityList;
    gridRows:();
    i0:0; while[i0<count ivL;
    i1:0; while[i1<count lvL;
    i2:0; while[i2<count mrL;
    i3:0; while[i3<count vvL;
    i4:0; while[i4<count coL;
    i5:0; while[i5<count jiL;
    i6:0; while[i6<count jmL;
    i7:0; while[i7<count jvL;
        gridRows:gridRows,enlist `initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility!(
            ivL i0;lvL i1;mrL i2;vvL i3;coL i4;jiL i5;jmL i6;jvL i7);
    i7+:1]; i6+:1]; i5+:1]; i4+:1]; i3+:1]; i2+:1]; i1+:1]; i0+:1];
    gridRows
 };

.bates.calibrateGrid:{[optionTable;marketData;paramGridDict;configDict]
    .calibration.validateMarketOptionTable optionTable;
    paramGrid:.bates.parameterGrid paramGridDict;
    mcConfig:$[`mcConfig in key configDict;configDict`mcConfig;.montecarlo.defaultMcConfig[]];
    calibrationTable:();
    gridIdx:0;
    while[gridIdx<count paramGrid;
        paramRow:paramGrid gridIdx;
        fullParams:paramRow,`riskFreeRate`dividendYield!(marketData`riskFreeRate;marketData`dividendYield);
        modelPrices:();
        marketPrices:();
        failedCount:0;
        optIdx:0;
        while[optIdx<count optionTable;
            optRow:optionTable optIdx;
            trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
                optRow`optionId;optRow`underlying;`equityOption;`european;optRow`optionType;optRow`strike;optRow`expiry;1f);
            batesConfig:`batesParams`mcConfig!(fullParams;mcConfig);
            priceFn:{.bates.priceEuropean[x 0;x 1;x 2]};
            priceResult:@[priceFn;(trade;marketData;batesConfig);{x}];
            if[10h=type priceResult;
                failedCount+:1];
            if[not 10h=type priceResult;
                modelPrices:modelPrices,priceResult`unitPrice;
                marketPrices:marketPrices,optRow`marketPrice];
            optIdx+:1];
        pricedCount:count modelPrices;
        rmseVal:$[pricedCount>0;.objective.rmse[modelPrices;marketPrices];0Nf];
        maeVal:$[pricedCount>0;.objective.mae[modelPrices;marketPrices];0Nf];
        diffs:modelPrices-marketPrices;
        sseVal:$[pricedCount>0;sum diffs*diffs;0Nf];
        calibrationTable:calibrationTable,enlist `calibrationId`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`rmse`mae`totalSse`pricedRows`failedRows`status`errorMessage!(
            gridIdx+1;paramRow`initialVariance;paramRow`longRunVariance;paramRow`meanReversion;paramRow`volOfVol;paramRow`correlation;
            paramRow`jumpIntensity;paramRow`jumpMean;paramRow`jumpVolatility;
            rmseVal;maeVal;sseVal;pricedCount;failedCount;`OK;"");
        gridIdx+:1];
    calibrationTable
 };

.bates.bestCalibration:{[calibrationTable]
    statusCol:calibrationTable`status;
    okMask:statusCol=`OK;
    okRows:calibrationTable where okMask;
    if[0=count okRows; '"No successful Bates calibration results"];
    rmseValues:okRows`rmse;
    bestIdx:rmseValues?min rmseValues;
    okRows bestIdx
 };
