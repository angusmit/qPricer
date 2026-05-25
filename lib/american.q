/ american.q - American put analysis: early exercise boundary extraction
/ The exercise boundary is the highest spot where immediate exercise is optimal.
/ Estimated from the FDM grid by comparing option value with intrinsic value.

.american.__exerciseTolerance:0.01;

/ --- Public ---

.american.extractEarlyExerciseBoundary:{[trade;marketData;model;config]
    .american.__validateAmericanPutTrade trade;
    .market.validateFlatMarketData marketData;
    .model.validateModel model;
    .config.validateFiniteDifferenceConfig config;
    gridResult:.engine.priceOptionWithGrid[trade;marketData;model;config];
    solverResult:gridResult`solverResult;
    spotGrid:solverResult`spotGrid;
    timeGrid:solverResult`timeGrid;
    valueGrid:solverResult`valueGrid;
    strike:trade`strike;
    expiry:trade`expiry;
    intrinsicValues:0f|strike-spotGrid;
    boundaryRows:.american.__extractBoundaryAtTime[spotGrid;timeGrid;valueGrid;intrinsicValues;trade;] each til count timeGrid;
    .american.__buildBoundaryTable boundaryRows
 };

.american.analyzeAmericanPut:{[trade;marketData;model;config]
    .american.__validateAmericanPutTrade trade;
    americanPriceResult:.engine.priceOption[trade;marketData;model;config];
    americanUnitPrice:americanPriceResult`unitPrice;
    europeanTrade:@[trade;`exerciseStyle;:;`european];
    europeanPriceResult:.engine.priceOption[europeanTrade;marketData;model;config];
    europeanUnitPrice:europeanPriceResult`unitPrice;
    premiumResult:`americanUnitPrice`europeanUnitPrice`earlyExercisePremium!(
        americanUnitPrice;europeanUnitPrice;americanUnitPrice-europeanUnitPrice);
    boundaryTable:.american.extractEarlyExerciseBoundary[trade;marketData;model;config];
    `priceResult`earlyExercisePremium`earlyExerciseBoundary!(
        americanPriceResult;premiumResult;boundaryTable)
 };

/ --- Internal ---

.american.__validateAmericanPutTrade:{[trade]
    .product.validateOptionTrade trade;
    if[not trade[`exerciseStyle]~`american;
        '"Early exercise boundary requires exerciseStyle = `american, got ",string trade`exerciseStyle];
    if[not trade[`optionType]~`put;
        '"Early exercise boundary is currently only supported for American put options"];
 };

.american.__extractBoundaryAtTime:{[spotGrid;timeGrid;valueGrid;intrinsicValues;trade;timeIdx]
    optionValuesAtTime:valueGrid[;timeIdx];
    timePoint:timeGrid timeIdx;
    remainingTime:trade[`expiry]-timePoint;
    / Find spots where intrinsic > 0 and option value is close to intrinsic
    hasIntrinsic:intrinsicValues>0f;
    priceDiff:abs optionValuesAtTime-intrinsicValues;
    exercisingMask:hasIntrinsic & priceDiff<=.american.__exerciseTolerance;
    exerciseSpots:spotGrid where exercisingMask;
    hasBoundary:0<count exerciseSpots;
    boundarySpot:0Nf;
    if[hasBoundary; boundarySpot:max exerciseSpots];
    `tradeId`underlying`optionType`timePoint`remainingTime`exerciseBoundary`hasExerciseRegion!(
        trade`tradeId;trade`underlying;trade`optionType;
        timePoint;remainingTime;boundarySpot;hasBoundary)
 };

.american.__buildBoundaryTable:{[boundaryRows]
    columnNames:`tradeId`underlying`optionType`timePoint`remainingTime`exerciseBoundary`hasExerciseRegion;
    flip columnNames!flip value each boundaryRows
 };
