/ american.q - American option analysis: early exercise boundary extraction (v0.14)
/ Supports both American put and American call.
/ Put boundary: highest spot where exercise is optimal.
/ Call boundary: lowest spot where exercise is optimal.

.american.__exerciseTolerance:0.01;

/ --- Public ---

.american.extractEarlyExerciseBoundary:{[trade;marketData;model;config]
    .american.__validateAmericanTrade trade;
    .market.validateFlatMarketData marketData;
    .model.validateModel model;
    .config.validateFiniteDifferenceConfig config;
    gridResult:.engine.priceOptionWithGrid[trade;marketData;model;config];
    solverResult:gridResult`solverResult;
    spotGrid:solverResult`spotGrid;
    timeGrid:solverResult`timeGrid;
    valueGrid:solverResult`valueGrid;
    intrinsicValues:.payoff.calculateIntrinsicValue[trade;spotGrid];
    boundaryRows:.american.__extractBoundaryAtTime[spotGrid;timeGrid;valueGrid;intrinsicValues;trade;] each til count timeGrid;
    .american.__buildBoundaryTable boundaryRows
 };

.american.analyzeAmericanOption:{[trade;marketData;model;config]
    .american.__validateAmericanTrade trade;
    americanPriceResult:.engine.priceOption[trade;marketData;model;config];
    americanUnitPrice:americanPriceResult`unitPrice;
    europeanTrade:@[trade;`exerciseStyle;:;`european];
    europeanPriceResult:.engine.priceOption[europeanTrade;marketData;model;config];
    europeanUnitPrice:europeanPriceResult`unitPrice;
    premiumResult:`americanUnitPrice`europeanUnitPrice`earlyExercisePremium!(
        americanUnitPrice;europeanUnitPrice;americanUnitPrice-europeanUnitPrice);
    `priceResult`earlyExercisePremium!(americanPriceResult;premiumResult)
 };

/ Backward-compatible alias
.american.analyzeAmericanPut:{[trade;marketData;model;config]
    .american.__validateAmericanPutTrade trade;
    analysisResult:.american.analyzeAmericanOption[trade;marketData;model;config];
    boundaryTable:.american.extractEarlyExerciseBoundary[trade;marketData;model;config];
    `priceResult`earlyExercisePremium`earlyExerciseBoundary!(
        analysisResult`priceResult;analysisResult`earlyExercisePremium;boundaryTable)
 };

/ --- Internal ---

.american.__validateAmericanTrade:{[trade]
    .product.validateOptionTrade trade;
    if[not trade[`exerciseStyle]~`american;
        '"Early exercise analysis requires exerciseStyle = `american, got ",string trade`exerciseStyle];
 };

.american.__validateAmericanPutTrade:{[trade]
    .american.__validateAmericanTrade trade;
    if[not trade[`optionType]~`put;
        '"analyzeAmericanPut requires optionType = `put"];
 };

.american.__extractBoundaryAtTime:{[spotGrid;timeGrid;valueGrid;intrinsicValues;trade;timeIdx]
    optionValuesAtTime:valueGrid[;timeIdx];
    timePoint:timeGrid timeIdx;
    remainingTime:trade[`expiry]-timePoint;
    hasIntrinsic:intrinsicValues>0f;
    priceDiff:abs optionValuesAtTime-intrinsicValues;
    exercisingMask:hasIntrinsic & priceDiff<=.american.__exerciseTolerance;
    exerciseSpots:spotGrid where exercisingMask;
    hasBoundary:0<count exerciseSpots;
    boundarySpot:0Nf;
    isCallOption:trade[`optionType]~`call;
    / Call: lowest spot in exercise region; Put: highest spot
    if[hasBoundary; boundarySpot:$[isCallOption;min exerciseSpots;max exerciseSpots]];
    `tradeId`underlying`optionType`timePoint`remainingTime`exerciseBoundary`hasExerciseRegion!(
        trade`tradeId;trade`underlying;trade`optionType;
        timePoint;remainingTime;boundarySpot;hasBoundary)
 };

.american.__buildBoundaryTable:{[boundaryRows]
    columnNames:`tradeId`underlying`optionType`timePoint`remainingTime`exerciseBoundary`hasExerciseRegion;
    flip columnNames!flip value each boundaryRows
 };
