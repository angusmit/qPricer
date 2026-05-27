/ merton.q - Merton jump-diffusion model (v0.24)

/ --- Validation ---

.merton.validateParams:{[mertonParams]
    if[not mertonParams[`volatility]>0f; '"volatility must be positive"];
    if[mertonParams[`jumpIntensity]<0f; '"jumpIntensity must be non-negative"];
    if[mertonParams[`jumpVolatility]<0f; '"jumpVolatility must be non-negative"];
 };

/ --- Jump compensator ---

.merton.jumpCompensator:{[mertonParams]
    muJ:mertonParams`jumpMean;
    sigmaJ:mertonParams`jumpVolatility;
    halfSigmaJ2:0.5*sigmaJ*sigmaJ;
    (exp muJ+halfSigmaJ2)-1f
 };

/ --- Series pricing (Merton 1976) ---

.merton.priceEuropeanSeries:{[optionType;spotVal;strikeVal;expiryVal;mertonParams;termCountValue]
    .merton.validateParams mertonParams;
    sigmaVal:mertonParams`volatility;
    lambdaVal:mertonParams`jumpIntensity;
    muJVal:mertonParams`jumpMean;
    sigmaJVal:mertonParams`jumpVolatility;
    riskFreeRate:mertonParams`riskFreeRate;
    dividendYield:mertonParams`dividendYield;
    halfSigmaJ2:0.5*sigmaJVal*sigmaJVal;
    kappaJ:(exp muJVal+halfSigmaJ2)-1f;
    lambdaPrime:lambdaVal*exp muJVal+halfSigmaJ2;
    lambdaPrimeT:lambdaPrime*expiryVal;
    logOnePlusKappa:muJVal+halfSigmaJ2;
    sigma2:sigmaVal*sigmaVal;
    sigmaJ2:sigmaJVal*sigmaJVal;
    lambdaKappa:lambdaVal*kappaJ;
    totalPrice:0f;
    termIdx:0;
    while[termIdx<termCountValue;
        logNumerator:termIdx*log lambdaPrimeT|1e-300;
        logFactorial:sum log 1f+til termIdx;
        logWeight:(logNumerator-logFactorial)-lambdaPrimeT;
        poissonW:exp logWeight;
        nSigmaJ2:termIdx*sigmaJ2;
        adjVar:sigma2+nSigmaJ2%expiryVal;
        adjSigma:sqrt adjVar;
        adjRate:(riskFreeRate-lambdaKappa)+termIdx*logOnePlusKappa%expiryVal;
        bsPrice:.validation.blackScholesClosedForm[optionType;spotVal;strikeVal;expiryVal;adjRate;dividendYield;adjSigma];
        totalPrice+:poissonW*bsPrice;
        termIdx+:1];
    totalPrice
 };

/ --- Monte Carlo paths with jumps ---

.merton.__poissonFromUniforms:{[lambdaDt;uniformVector]
    prob:exp neg lambdaDt;
    cumProb:prob;
    nValues:(count uniformVector)#0;
    termIdx:1;
    while[(termIdx<20) and cumProb<0.9999999;
        nValues:nValues+uniformVector>=cumProb;
        prob*:lambdaDt%termIdx;
        cumProb+:prob;
        termIdx+:1];
    nValues
 };

.merton.simulateJumpDiffusionPaths:{[spotVal;mertonParams;expiry;mcConfig]
    .merton.validateParams mertonParams;
    pathCountValue:mcConfig`pathCount;
    timeStepCountValue:mcConfig`timeStepCount;
    randomSeedValue:mcConfig`randomSeed;
    sigmaVal:mertonParams`volatility;
    lambdaVal:mertonParams`jumpIntensity;
    muJVal:mertonParams`jumpMean;
    sigmaJVal:mertonParams`jumpVolatility;
    riskFreeRate:mertonParams`riskFreeRate;
    dividendYield:mertonParams`dividendYield;
    dtVal:expiry%timeStepCountValue;
    sqrtDt:sqrt dtVal;
    lambdaDt:lambdaVal*dtVal;
    halfSigma2:0.5*sigmaVal*sigmaVal;
    halfSigmaJ2:0.5*sigmaJVal*sigmaJVal;
    kappaJ:(exp muJVal+halfSigmaJ2)-1f;
    lambdaKappa:lambdaVal*kappaJ;
    rMinusQ:riskFreeRate-dividendYield;
    compensatedDrift:(rMinusQ-lambdaKappa)-halfSigma2;
    driftPerStep:compensatedDrift*dtVal;
    totalNormals:pathCountValue*timeStepCountValue;
    diffNormals:.montecarlo.__generateNormals[totalNormals;randomSeedValue];
    diffMatrix:timeStepCountValue cut diffNormals;
    jumpNormals:.montecarlo.__generateNormals[totalNormals;randomSeedValue+1];
    jumpMatrix:timeStepCountValue cut jumpNormals;
    system "S ",string randomSeedValue+2;
    allUniforms:totalNormals?1f;
    uniformMatrix:timeStepCountValue cut allUniforms;
    currentSpot:pathCountValue#spotVal;
    totalJumpCounts:pathCountValue#0;
    spotPathList:enlist currentSpot;
    stepIdx:0;
    while[stepIdx<timeStepCountValue;
        zDiff:diffMatrix[;stepIdx];
        zJump:jumpMatrix[;stepIdx];
        uPoisson:uniformMatrix[;stepIdx];
        jumpCounts:.merton.__poissonFromUniforms[lambdaDt;uPoisson];
        totalJumpCounts+:jumpCounts;
        floatJumps:`float$jumpCounts;
        sqrtJumps:sqrt floatJumps;
        jumpSums:(floatJumps*muJVal)+sqrtJumps*sigmaJVal*zJump;
        diffPart:driftPerStep+sigmaVal*sqrtDt*zDiff;
        logReturn:diffPart+jumpSums;
        currentSpot:currentSpot*exp logReturn;
        spotPathList:spotPathList,enlist currentSpot;
        stepIdx+:1];
    spotPaths:flip spotPathList;
    `pathMatrix`finalSpot`totalJumpCount`status`errorMessage!(
        spotPaths;currentSpot;totalJumpCounts;`OK;"")
 };

/ --- European pricing wrapper ---

.merton.priceEuropean:{[trade;marketData;configDict]
    mertonParams:configDict`mertonParams;
    pricingMethod:$[`pricingMethod in key configDict;configDict`pricingMethod;`series];
    termCountValue:$[`termCount in key configDict;configDict`termCount;30];
    spotVal:marketData`spot;
    expiryVal:trade`expiry;
    strikeVal:trade`strike;
    optionType:trade`optionType;
    if[pricingMethod~`series;
        seriesPrice:.merton.priceEuropeanSeries[optionType;spotVal;strikeVal;expiryVal;mertonParams;termCountValue];
        :`tradeId`underlying`optionType`productType`modelType`pricingMethod`unitPrice`notionalPrice`standardError`lowerConfidence`upperConfidence`status`statusMessage!(
            trade`tradeId;trade`underlying;optionType;`equityOption;`merton;`series;
            seriesPrice;seriesPrice*trade`notional;0Nf;0Nf;0Nf;`OK;"")];
    if[pricingMethod~`monteCarlo;
        mcConfig:$[`mcConfig in key configDict;configDict`mcConfig;.montecarlo.defaultMcConfig[]];
        pathResult:.merton.simulateJumpDiffusionPaths[spotVal;mertonParams;expiryVal;mcConfig];
        terminalSpots:pathResult`finalSpot;
        payoffVector:$[optionType~`call;0f|terminalSpots-strikeVal;0f|strikeVal-terminalSpots];
        priceResult:.montecarlo.priceFromPayoffs[payoffVector;mertonParams`riskFreeRate;expiryVal;mcConfig`confidenceLevel];
        :`tradeId`underlying`optionType`productType`modelType`pricingMethod`unitPrice`notionalPrice`standardError`lowerConfidence`upperConfidence`status`statusMessage!(
            trade`tradeId;trade`underlying;optionType;`equityOption;`merton;`monteCarlo;
            priceResult`price;priceResult[`price]*trade`notional;
            priceResult`standardError;priceResult`lowerConfidence;priceResult`upperConfidence;`OK;"")];
    '"Unsupported pricingMethod: ",string pricingMethod
 };

/ --- Greeks ---

.merton.bumpGreek:{[trade;marketData;configDict;greekName]
    baseResult:.merton.priceEuropean[trade;marketData;configDict];
    basePrice:baseResult`unitPrice;
    bumpedConfig:configDict;
    bumpedMkt:marketData;
    bumpSize:0f;
    mertonParams:configDict`mertonParams;
    if[greekName~`delta;
        bumpSize:marketData[`spot]*0.01;
        bumpedMkt:@[marketData;`spot;:;marketData[`spot]+bumpSize]];
    if[greekName~`vega;
        bumpSize:0.01;
        bumpedParams:@[mertonParams;`volatility;:;mertonParams[`volatility]+bumpSize];
        bumpedConfig:@[configDict;`mertonParams;:;bumpedParams]];
    if[greekName~`rho;
        bumpSize:0.0001;
        bumpedParams:@[mertonParams;`riskFreeRate;:;mertonParams[`riskFreeRate]+bumpSize];
        bumpedConfig:@[configDict;`mertonParams;:;bumpedParams]];
    if[greekName~`jumpIntensitySensitivity;
        bumpSize:0.01;
        bumpedParams:@[mertonParams;`jumpIntensity;:;mertonParams[`jumpIntensity]+bumpSize];
        bumpedConfig:@[configDict;`mertonParams;:;bumpedParams]];
    if[greekName~`jumpMeanSensitivity;
        bumpSize:0.01;
        bumpedParams:@[mertonParams;`jumpMean;:;mertonParams[`jumpMean]+bumpSize];
        bumpedConfig:@[configDict;`mertonParams;:;bumpedParams]];
    if[greekName~`jumpVolatilitySensitivity;
        bumpSize:0.01;
        bumpedParams:@[mertonParams;`jumpVolatility;:;mertonParams[`jumpVolatility]+bumpSize];
        bumpedConfig:@[configDict;`mertonParams;:;bumpedParams]];
    if[bumpSize=0f; '"Unsupported greekName: ",string greekName];
    bumpedResult:.merton.priceEuropean[trade;bumpedMkt;bumpedConfig];
    (bumpedResult[`unitPrice]-basePrice)%bumpSize
 };

/ --- Path diagnostics ---

.merton.pathDiagnostics:{[pathResult]
    finalSpotVector:pathResult`finalSpot;
    jumpCounts:pathResult`totalJumpCount;
    spotPaths:pathResult`pathMatrix;
    `pathCount`timeStepCount`averageFinalSpot`averageTotalJumpCount`maxTotalJumpCount`minFinalSpot`maxFinalSpot!(
        count finalSpotVector;
        (count spotPaths 0)-1;
        avg finalSpotVector;
        avg `float$jumpCounts;
        max jumpCounts;
        min finalSpotVector;
        max finalSpotVector)
 };

/ --- Calibration ---

.merton.parameterGrid:{[paramGridDict]
    volList:paramGridDict`volatilityList;
    liList:paramGridDict`jumpIntensityList;
    jmList:paramGridDict`jumpMeanList;
    jvList:paramGridDict`jumpVolatilityList;
    gridRows:();
    vIdx:0;
    while[vIdx<count volList;
        lIdx:0;
        while[lIdx<count liList;
            mIdx:0;
            while[mIdx<count jmList;
                jIdx:0;
                while[jIdx<count jvList;
                    gridRows:gridRows,enlist `volatility`jumpIntensity`jumpMean`jumpVolatility!(
                        volList vIdx;liList lIdx;jmList mIdx;jvList jIdx);
                    jIdx+:1];
                mIdx+:1];
            lIdx+:1];
        vIdx+:1];
    gridRows
 };

.merton.calibrateGrid:{[optionTable;marketData;paramGridDict;configDict]
    .calibration.validateMarketOptionTable optionTable;
    paramGrid:.merton.parameterGrid paramGridDict;
    calibrationTable:();
    gridIdx:0;
    while[gridIdx<count paramGrid;
        paramRow:paramGrid gridIdx;
        fullParams:paramRow,`riskFreeRate`dividendYield!(marketData`riskFreeRate;marketData`dividendYield);
        modelPrices:();
        marketPrices:();
        failedCount:0;
        termCountValue:$[`termCount in key configDict;configDict`termCount;30];
        optIdx:0;
        while[optIdx<count optionTable;
            optRow:optionTable optIdx;
            seriesPrice:@[.merton.priceEuropeanSeries[optRow`optionType;marketData`spot;optRow`strike;optRow`expiry;fullParams;];termCountValue;{x}];
            if[10h=type seriesPrice;
                failedCount+:1];
            if[not 10h=type seriesPrice;
                modelPrices:modelPrices,seriesPrice;
                marketPrices:marketPrices,optRow`marketPrice];
            optIdx+:1];
        pricedCount:count modelPrices;
        rmseVal:$[pricedCount>0;.objective.rmse[modelPrices;marketPrices];0Nf];
        maeVal:$[pricedCount>0;.objective.mae[modelPrices;marketPrices];0Nf];
        diffs:modelPrices-marketPrices;
        sseVal:$[pricedCount>0;sum diffs*diffs;0Nf];
        calibrationTable:calibrationTable,enlist `calibrationId`volatility`jumpIntensity`jumpMean`jumpVolatility`rmse`mae`totalSse`pricedRows`failedRows`status`errorMessage!(
            gridIdx+1;paramRow`volatility;paramRow`jumpIntensity;paramRow`jumpMean;paramRow`jumpVolatility;
            rmseVal;maeVal;sseVal;pricedCount;failedCount;`OK;"");
        gridIdx+:1];
    calibrationTable
 };

.merton.bestCalibration:{[calibrationTable]
    statusCol:calibrationTable`status;
    okMask:statusCol=`OK;
    okRows:calibrationTable where okMask;
    if[0=count okRows; '"No successful Merton calibration results"];
    rmseValues:okRows`rmse;
    bestIdx:rmseValues?min rmseValues;
    okRows bestIdx
 };
