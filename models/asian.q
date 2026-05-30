/ asian.q - Asian option pricing via Monte Carlo (v0.16)

/ --- Public ---

.asian.priceAsianOption:{[trade;marketData;configDict]
    .asian.validateAsianTrade trade;
    spotVal:marketData`spot;
    riskFreeRate:marketData`riskFreeRate;
    dividendYield:marketData`dividendYield;
    volatility:marketData`volatility;
    expiry:trade`expiry;
    observationCountValue:trade`observationCount;
    mcConfig:$[`mcConfig in key configDict; configDict`mcConfig; .montecarlo.defaultMcConfig[]];
    .montecarlo.validateMcConfig mcConfig;
    / Set timeStepCount = observationCount for Asian
    mcConfigAsian:@[mcConfig;`timeStepCount;:;observationCountValue];
    pathMatrix:.montecarlo.simulateGBMPaths[spotVal;riskFreeRate;dividendYield;volatility;expiry;mcConfigAsian];
    averageVector:.asian.averagePath[pathMatrix;trade`averageType];
    payoffVector:.asian.payoff[averageVector;trade`strike;trade`optionType];
    priceResult:.montecarlo.priceFromPayoffs[payoffVector;riskFreeRate;expiry;mcConfig`confidenceLevel];
    `tradeId`underlying`optionType`productType`unitPrice`notionalPrice`method`modelName`standardError`lowerConfidence`upperConfidence`status`statusMessage!(
        trade`tradeId;trade`underlying;trade`optionType;`asianOption;
        priceResult`price;priceResult[`price]*trade`notional;
        `monteCarlo;`asianOption;
        priceResult`standardError;priceResult`lowerConfidence;priceResult`upperConfidence;
        `OK;"")
 };

/ --- Path averaging ---

.asian.averagePath:{[pathMatrix;averageType]
    if[averageType~`arithmetic; :avg each pathMatrix];
    if[averageType~`geometric; :exp avg each log each pathMatrix];
    '"Unsupported averageType: ",string averageType
 };

.asian.payoff:{[averageVector;strike;optionType]
    if[optionType~`call; :0f|averageVector-strike];
    if[optionType~`put; :0f|strike-averageVector];
    '"Unsupported optionType: ",string optionType
 };

/ --- Geometric Asian closed-form (Kemna-Vorst) ---

.asian.geometricAsianClosedForm:{[optionType;spotVal;strike;expiry;riskFreeRate;dividendYield;volatility;observationCountValue]
    nVal:`float$observationCountValue;
    sigmaGSquared:volatility*volatility*expiry*(nVal+1f)*(2f*nVal+1f)%(6f*nVal*nVal);
    sigmaG:sqrt sigmaGSquared;
    halfVar:0.5*volatility*volatility;
    muAdjusted:((riskFreeRate-dividendYield)-halfVar)*(nVal+1f)%(2f*nVal)*expiry;
    forwardG:spotVal*exp muAdjusted+0.5*sigmaGSquared;
    d1Val:(log[forwardG%strike]+0.5*sigmaGSquared)%sigmaG;
    d2Val:d1Val-sigmaG;
    discountFactor:exp neg riskFreeRate*expiry;
    if[optionType~`call;
        :discountFactor*(forwardG*.montecarlo.__normalCDF[d1Val])-(strike*.montecarlo.__normalCDF d2Val)];
    if[optionType~`put;
        :discountFactor*(strike*.montecarlo.__normalCDF[neg d2Val])-(forwardG*.montecarlo.__normalCDF neg d1Val)];
    '"Unsupported optionType"
 };

/ --- Validation ---

.asian.validateAsianTrade:{[trade]
    if[not trade[`productType]~`asianOption; '"productType must be `asianOption"];
    if[not trade[`optionType] in `call`put; '"optionType must be `call or `put"];
    if[not trade[`exerciseStyle]~`european; '"Asian options only support European exercise"];
    if[not trade[`averageType] in `arithmetic`geometric; '"averageType must be `arithmetic or `geometric"];
    if[not trade[`strike]>0f; '"strike must be positive"];
    if[not trade[`expiry]>0f; '"expiry must be positive"];
    if[not trade[`observationCount]>0; '"observationCount must be positive"];
 };

/ --- Control variate: arithmetic Asian with geometric Asian control (v0.19) ---

.asian.asianPayoffPair:{[pathMatrix;trade]
    arithmeticAvg:avg each pathMatrix;
    geometricAvg:exp avg each log each pathMatrix;
    arithmeticPayoff:.asian.payoff[arithmeticAvg;trade`strike;trade`optionType];
    geometricPayoff:.asian.payoff[geometricAvg;trade`strike;trade`optionType];
    `arithmeticPayoff`geometricPayoff!(arithmeticPayoff;geometricPayoff)
 };

.asian.priceAsianOptionWithControlVariate:{[trade;marketData;configDict]
    .asian.validateAsianTrade trade;
    if[not trade[`averageType]~`arithmetic; '"Control variate only for arithmetic Asian"];
    spotVal:marketData`spot;
    riskFreeRate:marketData`riskFreeRate;
    dividendYield:marketData`dividendYield;
    volatility:marketData`volatility;
    expiry:trade`expiry;
    observationCountValue:trade`observationCount;
    mcConfig:$[`mcConfig in key configDict; configDict`mcConfig; .montecarlo.defaultMcConfig[]];
    .montecarlo.validateMcConfig mcConfig;
    mcConfigAsian:@[mcConfig;`timeStepCount;:;observationCountValue];
    pathMatrix:.montecarlo.simulateGBMPaths[spotVal;riskFreeRate;dividendYield;volatility;expiry;mcConfigAsian];
    payoffPair:.asian.asianPayoffPair[pathMatrix;trade];
    / Geometric Asian closed-form (undiscounted expected payoff)
    geoClosedForm:.asian.geometricAsianClosedForm[trade`optionType;spotVal;trade`strike;expiry;riskFreeRate;dividendYield;volatility;observationCountValue];
    discountFactor:exp neg riskFreeRate*expiry;
    controlExpectedPayoff:geoClosedForm%discountFactor;
    / Apply control variate
    adjustedPayoff:.variance.controlVariateAdjust[payoffPair`arithmeticPayoff;payoffPair`geometricPayoff;controlExpectedPayoff];
    betaValue:.variance.controlVariateBeta[payoffPair`arithmeticPayoff;payoffPair`geometricPayoff];
    / Price from adjusted payoffs
    priceResult:.montecarlo.priceFromPayoffs[adjustedPayoff;riskFreeRate;expiry;mcConfig`confidenceLevel];
    `tradeId`underlying`optionType`productType`unitPrice`notionalPrice`method`modelName`standardError`lowerConfidence`upperConfidence`betaValue`controlVariate`status`statusMessage!(
        trade`tradeId;trade`underlying;trade`optionType;`asianOption;
        priceResult`price;priceResult[`price]*trade`notional;
        `monteCarlo;`asianOptionCV;
        priceResult`standardError;priceResult`lowerConfidence;priceResult`upperConfidence;
        betaValue;`geometricAsian;`OK;"")
 };
