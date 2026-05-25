/ engine.q — pricing engine

/ --- Public ---

.engine.priceOption:{[trade;marketData;model;config]
    .engine.__validateInputs[trade;marketData;model;config];
    solverResult:.engine.__runSolver[trade;marketData;model;config];
    unitPrice:.engine.__interpolatePriceFromGrid[
        marketData`spot; solverResult`spotGrid; solverResult`valueAtTimeZero; config`interpolationMethod];
    notionalPrice:unitPrice * trade`notional;
    `tradeId`underlying`optionType`unitPrice`notionalPrice`method`modelName!(
        trade`tradeId; trade`underlying; trade`optionType;
        unitPrice; notionalPrice; config`method; model`modelName)
 };

.engine.priceOptionWithGrid:{[trade;marketData;model;config]
    configFull:@[config;`returnFullGrid;:;1b];
    .engine.__validateInputs[trade;marketData;model;configFull];
    solverResult:.engine.__runSolver[trade;marketData;model;configFull];
    unitPrice:.engine.__interpolatePriceFromGrid[
        marketData`spot; solverResult`spotGrid; solverResult`valueAtTimeZero; configFull`interpolationMethod];
    notionalPrice:unitPrice * trade`notional;
    priceResult:`tradeId`underlying`optionType`unitPrice`notionalPrice`method`modelName!(
        trade`tradeId; trade`underlying; trade`optionType;
        unitPrice; notionalPrice; configFull`method; model`modelName);
    `priceResult`solverResult!(priceResult;solverResult)
 };

/ --- Internal ---

.engine.__validateInputs:{[trade;marketData;model;config]
    .product.validateOptionTrade[trade];
    .market.validateFlatMarketData[marketData];
    .model.validateModel[model];
    .config.validateFiniteDifferenceConfig[config];
    if[not trade[`underlying] ~ marketData`underlying;
        '"Underlying mismatch: trade has ",string[trade`underlying]," but marketData has ",string marketData`underlying];
 };

.engine.__runSolver:{[trade;marketData;model;config]
    solver:.engine.__selectSolver[config];
    solver[trade;marketData;model;config]
 };

.engine.__selectSolver:{[config]
    if[config[`method]~`explicit; :.solver.solveExplicitFiniteDifference];
    '"Unsupported solver method: ",string config`method
 };

.engine.__interpolatePriceFromGrid:{[spot;spotGrid;valueAtTimeZero;interpolationMethod]
    if[interpolationMethod~`linear; :.utilities.linearInterpolate[spot;spotGrid;valueAtTimeZero]];
    '"Unsupported interpolation method: ",string interpolationMethod
 };
