/ lookback.q - lookback option pricing via Monte Carlo (v0.18)

/ --- Public ---

.lookback.priceLookbackOption:{[trade;marketData;configDict]
    .lookback.validateLookbackTrade trade;
    spotVal:marketData`spot;
    riskFreeRate:marketData`riskFreeRate;
    dividendYield:marketData`dividendYield;
    volatility:marketData`volatility;
    expiry:trade`expiry;
    mcConfig:$[`mcConfig in key configDict; configDict`mcConfig; .montecarlo.defaultMcConfig[]];
    .montecarlo.validateMcConfig mcConfig;
    mcConfigLB:@[mcConfig;`timeStepCount;:;trade`observationCount];
    pathMatrix:.montecarlo.simulateGBMPaths[spotVal;riskFreeRate;dividendYield;volatility;expiry;mcConfigLB];
    payoffVector:.lookback.pathPayoff[pathMatrix;trade];
    priceResult:.montecarlo.priceFromPayoffs[payoffVector;riskFreeRate;expiry;mcConfig`confidenceLevel];
    `tradeId`underlying`optionType`productType`unitPrice`notionalPrice`method`modelName`standardError`lowerConfidence`upperConfidence`status`statusMessage!(
        trade`tradeId;trade`underlying;trade`optionType;`lookbackOption;
        priceResult`price;priceResult[`price]*trade`notional;
        `monteCarlo;`lookbackOption;
        priceResult`standardError;priceResult`lowerConfidence;priceResult`upperConfidence;
        `OK;"")
 };

.lookback.pathPayoff:{[pathMatrix;trade]
    lookbackStyle:trade`lookbackStyle;
    optionType:trade`optionType;
    pathMaxVector:.pathdiag.pathMaximum pathMatrix;
    pathMinVector:.pathdiag.pathMinimum pathMatrix;
    pathFinalVector:.pathdiag.pathFinal pathMatrix;
    if[lookbackStyle~`fixed;
        strike:trade`strike;
        if[optionType~`call; :0f|pathMaxVector-strike];
        if[optionType~`put; :0f|strike-pathMinVector];
        '"Unsupported optionType for fixed lookback"];
    if[lookbackStyle~`floating;
        if[optionType~`call; :0f|pathFinalVector-pathMinVector];
        if[optionType~`put; :0f|pathMaxVector-pathFinalVector];
        '"Unsupported optionType for floating lookback"];
    '"Unsupported lookbackStyle: ",string lookbackStyle
 };

/ --- Greeks via bump-and-reprice ---

.lookback.bumpGreek:{[trade;marketData;configDict;greekName]
    baseResult:.lookback.priceLookbackOption[trade;marketData;configDict];
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
    bumpedResult:.lookback.priceLookbackOption[trade;bumpedMkt;configDict];
    (bumpedResult[`unitPrice]-basePrice)%bumpSize
 };

/ --- Validation ---

.lookback.validateLookbackTrade:{[trade]
    if[not trade[`productType]~`lookbackOption; '"productType must be `lookbackOption"];
    if[not trade[`lookbackStyle] in `fixed`floating; '"lookbackStyle must be `fixed or `floating"];
    if[not trade[`optionType] in `call`put; '"optionType must be `call or `put"];
    if[not trade[`exerciseStyle]~`european; '"Lookback options only support European exercise"];
    if[not trade[`expiry]>0f; '"expiry must be positive"];
    if[not trade[`observationCount]>0; '"observationCount must be positive"];
    if[trade[`lookbackStyle]~`fixed;
        if[not trade[`strike]>0f; '"fixed-strike lookback requires positive strike"]];
 };
