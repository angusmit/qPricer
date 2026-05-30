/ montecarlo.q - Monte Carlo simulation engine (v0.16)

/ --- Normal distribution utilities ---

.montecarlo.__normalCDF:{[xVal]
    if[xVal<-10f; :0f];
    if[xVal>10f; :1f];
    isNeg:xVal<0f;
    absX:abs xVal;
    kVal:1f%(1f+0.2316419*absX);
    nPdf:(1f%sqrt 6.28318530717959)*exp neg 0.5*absX*absX;
    polyVal:kVal*(0.319381530+kVal*(-0.356563782+kVal*(1.781477937+kVal*(-1.821255978+1.330274429*kVal))));
    resultVal:1f-nPdf*polyVal;
    $[isNeg;1f-resultVal;resultVal]
 };

.montecarlo.__zScoreForConfidence:{[confidenceLevel]
    if[(abs confidenceLevel-0.90)<0.001; :1.6449];
    if[(abs confidenceLevel-0.95)<0.001; :1.9600];
    if[(abs confidenceLevel-0.99)<0.001; :2.5758];
    1.9600
 };

/ --- Random number generation ---

.montecarlo.__generateNormals:{[totalNormals;randomSeedValue]
    system "S ",string randomSeedValue;
    pairCount:ceiling totalNormals%2;
    u1:0.0000001|pairCount?1f;
    u2:pairCount?1f;
    radius:sqrt -2f*log u1;
    angle:6.28318530717959*u2;
    normals:raze (radius*cos angle;radius*sin angle);
    totalNormals#normals
 };

.montecarlo.generateNormalPaths:{[pathCountValue;timeStepCountValue;randomSeedValue]
    totalNormals:pathCountValue*timeStepCountValue;
    flatNormals:.montecarlo.__generateNormals[totalNormals;randomSeedValue];
    timeStepCountValue cut flatNormals
 };

.montecarlo.generateAntitheticNormals:{[pathCountValue;timeStepCountValue;randomSeedValue]
    halfPaths:ceiling pathCountValue%2;
    baseMatrix:.montecarlo.generateNormalPaths[halfPaths;timeStepCountValue;randomSeedValue];
    antiMatrix:neg each baseMatrix;
    pathCountValue#(baseMatrix,antiMatrix)
 };

/ --- GBM simulation ---

.montecarlo.simulateGBMPaths:{[spotVal;riskFreeRate;dividendYield;volatility;expiry;mcConfig]
    pathCountValue:mcConfig`pathCount;
    timeStepCountValue:mcConfig`timeStepCount;
    randomSeedValue:mcConfig`randomSeed;
    useAntithetic:mcConfig`antithetic;
    dtVal:expiry%timeStepCountValue;
    halfVariance:0.5*volatility*volatility;
    driftPerStep:((riskFreeRate-dividendYield)-halfVariance)*dtVal;
    diffusionPerStep:volatility*sqrt dtVal;
    normalMatrix:.montecarlo.generateNormalPaths[pathCountValue;timeStepCountValue;randomSeedValue];
    if[useAntithetic; normalMatrix:.montecarlo.generateAntitheticNormals[pathCountValue;timeStepCountValue;randomSeedValue]];
    logIncrements:driftPerStep + diffusionPerStep * normalMatrix;
    logSpotPaths:(log spotVal) + sums each logIncrements;
    exp each logSpotPaths
 };

/ --- Pricing from payoffs ---

.montecarlo.discountedPayoff:{[payoffVector;riskFreeRate;expiry]
    (exp neg riskFreeRate*expiry)*payoffVector
 };

.montecarlo.priceFromPayoffs:{[payoffVector;riskFreeRate;expiry;confidenceLevel]
    discountedPayoffs:.montecarlo.discountedPayoff[payoffVector;riskFreeRate;expiry];
    priceVal:avg discountedPayoffs;
    pathCountValue:count payoffVector;
    stdDevVal:dev discountedPayoffs;
    standardError:stdDevVal%sqrt `float$pathCountValue;
    zScore:.montecarlo.__zScoreForConfidence confidenceLevel;
    `price`standardError`lowerConfidence`upperConfidence`pathCount!(
        priceVal;standardError;priceVal-zScore*standardError;priceVal+zScore*standardError;pathCountValue)
 };

/ --- MC European (for validation) ---

.montecarlo.priceEuropeanMC:{[optionType;spotVal;strike;expiry;riskFreeRate;dividendYield;volatility;mcConfig]
    pathMatrix:.montecarlo.simulateGBMPaths[spotVal;riskFreeRate;dividendYield;volatility;expiry;mcConfig];
    terminalSpots:last each pathMatrix;
    payoffVector:$[optionType~`call;0f|terminalSpots-strike;0f|strike-terminalSpots];
    .montecarlo.priceFromPayoffs[payoffVector;riskFreeRate;expiry;mcConfig`confidenceLevel]
 };

/ --- Bump-and-reprice Greeks ---

.montecarlo.bumpGreek:{[trade;marketData;configDict;greekName]
    baseResult:.asian.priceAsianOption[trade;marketData;configDict];
    basePrice:baseResult`unitPrice;
    bumpedMkt:marketData;
    bumpSize:0f;
    if[greekName~`delta;
        bumpSize:marketData[`spot]*0.01;
        bumpedMkt:@[marketData;`spot;:;marketData[`spot]+bumpSize]];
    if[greekName~`vega;
        bumpSize:0.01;
        bumpedMkt:@[marketData;`volatility;:;marketData[`volatility]+bumpSize]];
    if[greekName~`rho;
        bumpSize:0.0001;
        bumpedMkt:@[marketData;`riskFreeRate;:;marketData[`riskFreeRate]+bumpSize]];
    if[bumpSize=0f; '"Unsupported greekName: ",string greekName];
    bumpedResult:.asian.priceAsianOption[trade;bumpedMkt;configDict];
    (bumpedResult[`unitPrice]-basePrice)%bumpSize
 };

/ --- Config ---

.montecarlo.defaultMcConfig:{[]
    .cfg.mc
 };

.montecarlo.validateMcConfig:{[mcConfig]
    if[not mcConfig[`pathCount]>0; '"MC pathCount must be positive"];
    if[not mcConfig[`timeStepCount]>0; '"MC timeStepCount must be positive"];
    if[not mcConfig[`confidenceLevel]>0f; '"MC confidenceLevel must be positive"];
    if[not mcConfig[`confidenceLevel]<1f; '"MC confidenceLevel must be less than 1"];
 };

/ --- Correlated Monte Carlo (v0.17) ---

.montecarlo.generateCorrelatedNormals:{[pathCountValue;timeStepCountValue;correlationMatrix;randomSeedValue]
    numSymbols:count correlationMatrix;
    totalRows:pathCountValue*timeStepCountValue;
    / Generate independent normals: totalRows x numSymbols
    allNormals:.montecarlo.__generateNormals[totalRows*numSymbols;randomSeedValue];
    independentMatrix:numSymbols cut allNormals;  / totalRows x numSymbols
    / Apply Cholesky: Z = W * L^T
    choleskyL:.correlation.__cholesky correlationMatrix;
    independentMatrix mmu flip choleskyL
 };

.montecarlo.simulateCorrelatedGBMPaths:{[spotVector;riskFreeRateVector;dividendYieldVector;volatilityVector;expiry;correlationMatrix;mcConfig]
    pathCountValue:mcConfig`pathCount;
    timeStepCountValue:mcConfig`timeStepCount;
    randomSeedValue:mcConfig`randomSeed;
    numSymbols:count spotVector;
    dtVal:expiry%timeStepCountValue;
    / Correlated normals: (pathCount*timeStepCount) x numSymbols
    correlatedNormals:.montecarlo.generateCorrelatedNormals[pathCountValue;timeStepCountValue;correlationMatrix;randomSeedValue];
    / Build per-symbol drift and diffusion
    halfVariances:0.5*volatilityVector*volatilityVector;
    driftPerStepVector:((riskFreeRateVector-dividendYieldVector)-halfVariances)*dtVal;
    diffusionPerStepVector:volatilityVector*sqrt dtVal;
    / For each symbol, extract its column of correlated normals and build paths
    pathData:numSymbols#enlist();
    symIdx:0;
    while[symIdx<numSymbols;
        symbolNormals:correlatedNormals[;symIdx];  / totalRows column
        symbolNormalMatrix:timeStepCountValue cut symbolNormals;  / pathCount x timeStepCount
        logIncrements:driftPerStepVector[symIdx]+diffusionPerStepVector[symIdx]*symbolNormalMatrix;
        logSpotPaths:(log spotVector symIdx)+sums each logIncrements;
        pathData[symIdx]:exp each logSpotPaths;
        symIdx+:1];
    pathData  / list of numSymbols pathMatrices, each pathCount x timeStepCount
 };
