/ engine.q - pricing engine (v0.14: knock-in via parity)

.engine.priceOption:{[trade;marketData;model;config]
    .engine.__validateInputs[trade;marketData;model;config];
    if[.product.isKnockIn trade;
        :.engine.__priceKnockInViaParity[trade;marketData;model;config]];
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
    if[.product.isKnockIn trade;
        :.engine.__priceKnockInWithGrid[trade;marketData;model;configFull]];
    solverResult:.engine.__runSolver[trade;marketData;model;configFull];
    unitPrice:.engine.__interpolatePriceFromGrid[
        marketData`spot; solverResult`spotGrid; solverResult`valueAtTimeZero; configFull`interpolationMethod];
    notionalPrice:unitPrice * trade`notional;
    priceResult:`tradeId`underlying`optionType`unitPrice`notionalPrice`method`modelName!(
        trade`tradeId; trade`underlying; trade`optionType;
        unitPrice; notionalPrice; configFull`method; model`modelName);
    `priceResult`solverResult!(priceResult;solverResult)
 };

/ --- Knock-in via parity: knockIn = vanilla - knockOut ---

.engine.__priceKnockInViaParity:{[trade;marketData;model;config]
    vanillaTrade:@[trade;`barrierType;:;`none];
    knockOutType:.product.knockInToKnockOut .product.getBarrierType trade;
    knockOutTrade:@[trade;`barrierType;:;knockOutType];
    vanillaResult:.engine.priceOption[vanillaTrade;marketData;model;config];
    knockOutResult:.engine.priceOption[knockOutTrade;marketData;model;config];
    knockInUnitPrice:vanillaResult[`unitPrice]-knockOutResult`unitPrice;
    knockInNotionalPrice:knockInUnitPrice*trade`notional;
    `tradeId`underlying`optionType`unitPrice`notionalPrice`method`modelName!(
        trade`tradeId;trade`underlying;trade`optionType;
        knockInUnitPrice;knockInNotionalPrice;config`method;model`modelName)
 };

.engine.__priceKnockInWithGrid:{[trade;marketData;model;config]
    vanillaTrade:@[trade;`barrierType;:;`none];
    knockOutType:.product.knockInToKnockOut .product.getBarrierType trade;
    knockOutTrade:@[trade;`barrierType;:;knockOutType];
    vanillaGrid:.engine.__runSolver[vanillaTrade;marketData;model;config];
    knockOutGrid:.engine.__runSolver[knockOutTrade;marketData;model;config];
    knockInValues:vanillaGrid[`valueAtTimeZero]-knockOutGrid`valueAtTimeZero;
    knockInValueGrid:vanillaGrid[`valueGrid]-knockOutGrid`valueGrid;
    knockInUnitPrice:.engine.__interpolatePriceFromGrid[
        marketData`spot;vanillaGrid`spotGrid;knockInValues;config`interpolationMethod];
    knockInNotionalPrice:knockInUnitPrice*trade`notional;
    priceResult:`tradeId`underlying`optionType`unitPrice`notionalPrice`method`modelName!(
        trade`tradeId;trade`underlying;trade`optionType;
        knockInUnitPrice;knockInNotionalPrice;config`method;model`modelName);
    solverResult:`spotGrid`timeGrid`valueGrid`valueAtTimeZero`metadata!(
        vanillaGrid`spotGrid;vanillaGrid`timeGrid;knockInValueGrid;knockInValues;vanillaGrid`metadata);
    `priceResult`solverResult!(priceResult;solverResult)
 };

/ --- Internal ---

.engine.__validateInputs:{[trade;marketData;model;config]
    .product.validateOptionTrade[trade];
    .market.validateMarketData[marketData;model];
    .model.validateModel[model];
    .config.validateFiniteDifferenceConfig[config];
    if[not trade[`underlying] ~ marketData`underlying;
        '"Underlying mismatch: trade has ",string[trade`underlying]," but marketData has ",string marketData`underlying];
 };

.engine.__runSolver:{[trade;marketData;model;config]
    solver:.engine.__selectSolver config;
    solver[trade;marketData;model;config]
 };

.engine.__selectSolver:{[config]
    if[config[`method]~`explicit; :.solver.solveExplicitFiniteDifference];
    if[config[`method]~`crankNicolson; :.solver.solveCrankNicolson];
    '"Unsupported solver method: ",string config`method
 };

.engine.__interpolatePriceFromGrid:{[spot;spotGrid;valueAtTimeZero;interpolationMethod]
    if[interpolationMethod~`linear; :.utilities.linearInterpolate[spot;spotGrid;valueAtTimeZero]];
    '"Unsupported interpolation method: ",string interpolationMethod
 };
