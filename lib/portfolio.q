/ portfolio.q - table-based portfolio pricing and batch risk
/ enlist dict in q creates a 1-row table; successive , builds the result.

.portfolio.__requiredColumns:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate;

/ --- Public ---

.portfolio.priceTrades:{[tradeTable;marketData;model;config]
    .portfolio.__validateTradeTable tradeTable;
    numTrades:count tradeTable;
    resultTable:();
    loopIdx:0;
    while[loopIdx<numTrades;
        singleResult:.portfolio.__priceSingleTradeSafely[tradeTable loopIdx;marketData;model;config];
        resultTable:resultTable,enlist singleResult;
        loopIdx+:1];
    resultTable
 };

.portfolio.calculatePortfolioGreeks:{[tradeTable;marketData;model;config]
    .portfolio.__validateTradeTable tradeTable;
    numTrades:count tradeTable;
    resultTable:();
    loopIdx:0;
    while[loopIdx<numTrades;
        singleResult:.portfolio.__calculateSingleTradeGreeksSafely[tradeTable loopIdx;marketData;model;config];
        resultTable:resultTable,enlist singleResult;
        loopIdx+:1];
    resultTable
 };

.portfolio.generatePortfolioScenarioReport:{[tradeTable;marketData;model;config]
    .portfolio.__validateTradeTable tradeTable;
    numTrades:count tradeTable;
    resultTable:();
    loopIdx:0;
    while[loopIdx<numTrades;
        singleBlock:.portfolio.__generateSingleTradeScenarioSafely[tradeTable loopIdx;marketData;model;config];
        resultTable:resultTable,singleBlock;
        loopIdx+:1];
    resultTable
 };

.portfolio.summarizePortfolioRisk:{[portfolioScenarioTable]
    okRows:portfolioScenarioTable where portfolioScenarioTable[`status]=`OK;
    scenarioNames:distinct okRows`scenario;
    resultTable:();
    scenIdx:0;
    while[scenIdx<count scenarioNames;
        scenName:scenarioNames scenIdx;
        scRows:okRows where okRows[`scenario]=scenName;
        summaryRow:`scenario`totalNotionalPrice`totalNotionalPnL!(scenName;sum scRows`notionalPrice;sum scRows`notionalPnL);
        resultTable:resultTable,enlist summaryRow;
        scenIdx+:1];
    resultTable
 };

/ --- Internal ---

.portfolio.__validateTradeTable:{[tradeTable]
    if[0=count tradeTable; '"Portfolio trade table is empty"];
    tableColumnNames:cols tradeTable;
    missingColumns:.portfolio.__requiredColumns where not .portfolio.__requiredColumns in tableColumnNames;
    if[0<count missingColumns;
        '"Missing required portfolio columns: ",", " sv string missingColumns];
 };

.portfolio.__priceSingleTradeSafely:{[tradeDictionary;marketData;model;config]
    barrierType:.product.getBarrierType tradeDictionary;
    / Route Asian to MC engine
    if[tradeDictionary[`productType]~`asianOption;
        :.portfolio.__priceAsianSafely[tradeDictionary;marketData;model;config]];
    / Route Basket to correlated MC engine
    if[tradeDictionary[`productType]~`basketOption;
        :.portfolio.__priceBasketSafely[tradeDictionary;marketData;model;config]];
    / Route Lookback to MC engine
    if[tradeDictionary[`productType]~`lookbackOption;
        :.portfolio.__priceLookbackSafely[tradeDictionary;marketData;model;config]];
    / Route Heston equityOption to Heston MC
    if[`modelType in key config;
        if[config[`modelType]~`heston;
            if[tradeDictionary[`productType]~`equityOption;
                :.portfolio.__priceHestonSafely[tradeDictionary;marketData;model;config]]];
        if[config[`modelType]~`merton;
            if[tradeDictionary[`productType]~`equityOption;
                :.portfolio.__priceMertonSafely[tradeDictionary;marketData;model;config]]]];
    functionResult:.[.engine.priceOption;(tradeDictionary;marketData;model;config);{x}];
    if[10h=type functionResult;
        :`tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            0Nf;0Nf;config`method;model`modelName;`ERROR;functionResult)];
    `tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
        tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
        tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
        functionResult`unitPrice;functionResult`notionalPrice;
        functionResult`method;functionResult`modelName;`OK;"")
 };

.portfolio.__priceAsianSafely:{[tradeDictionary;marketData;model;config]
    barrierType:`none;
    mcConfig:$[`mcConfig in key config; config`mcConfig; .montecarlo.defaultMcConfig[]];
    asianConfig:`mcConfig`model`fdmConfig!(mcConfig;model;config);
    functionResult:.[.asian.priceAsianOption;(tradeDictionary;marketData;asianConfig);{x}];
    if[10h=type functionResult;
        :`tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            0Nf;0Nf;`monteCarlo;`asianOption;`ERROR;functionResult)];
    `tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
        tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
        tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
        functionResult`unitPrice;functionResult`notionalPrice;
        `monteCarlo;`asianOption;`OK;"")
 };

.portfolio.__priceBasketSafely:{[tradeDictionary;marketData;model;config]
    barrierType:`none;
    / Basket needs marketDataBook and correlationTable from config
    / If marketData is flat (not a book), wrap it
    correlationTable:$[`correlationTable in key config; config`correlationTable; '"No correlationTable in config for basket option"];
    / For single-symbol market data, basket can't route — needs market data book
    functionResult:.[{[td;md;cfg;ct]
        .basket.priceBasketOption[td;md;ct;cfg]
     };(tradeDictionary;marketData;config;correlationTable);{x}];
    if[10h=type functionResult;
        :`tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary[`basketSymbols]0;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            0Nf;0Nf;`monteCarlo;`basketOption;`ERROR;functionResult)];
    `tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
        tradeDictionary`tradeId;tradeDictionary[`basketSymbols]0;tradeDictionary`productType;
        tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
        functionResult`unitPrice;functionResult`notionalPrice;
        `monteCarlo;`basketOption;`OK;"")
 };

.portfolio.__priceLookbackSafely:{[tradeDictionary;marketData;model;config]
    barrierType:`none;
    mcConfig:$[`mcConfig in key config; config`mcConfig; .montecarlo.defaultMcConfig[]];
    lookbackConfig:`mcConfig`model`fdmConfig!(mcConfig;model;config);
    functionResult:.[.lookback.priceLookbackOption;(tradeDictionary;marketData;lookbackConfig);{x}];
    if[10h=type functionResult;
        :`tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            0Nf;0Nf;`monteCarlo;`lookbackOption;`ERROR;functionResult)];
    `tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
        tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
        tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
        functionResult`unitPrice;functionResult`notionalPrice;
        `monteCarlo;`lookbackOption;`OK;"")
 };

.portfolio.__priceHestonSafely:{[tradeDictionary;marketData;model;config]
    barrierType:`none;
    if[not tradeDictionary[`exerciseStyle]~`european;
        :`tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            0Nf;0Nf;`monteCarlo;`heston;`ERROR;"Heston only supports European exercise")];
    functionResult:.[.heston.priceEuropean;(tradeDictionary;marketData;config);{x}];
    if[10h=type functionResult;
        :`tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            0Nf;0Nf;`monteCarlo;`heston;`ERROR;functionResult)];
    `tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
        tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
        tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
        functionResult`unitPrice;functionResult`notionalPrice;
        `monteCarlo;`heston;`OK;"")
 };

.portfolio.__priceMertonSafely:{[tradeDictionary;marketData;model;config]
    barrierType:`none;
    if[not tradeDictionary[`exerciseStyle]~`european;
        :`tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            0Nf;0Nf;`series;`merton;`ERROR;"Merton only supports European exercise")];
    mertonFn:{.merton.priceEuropean[x 0;x 1;x 2]};
    functionResult:@[mertonFn;(tradeDictionary;marketData;config);{x}];
    if[10h=type functionResult;
        :`tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            0Nf;0Nf;`series;`merton;`ERROR;functionResult)];
    `tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
        tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
        tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
        functionResult`unitPrice;functionResult`notionalPrice;
        `series;`merton;`OK;"")
 };

.portfolio.__isGreeksSupported:{[tradeDictionary]
    / v0.14: Greeks via bump-and-reprice work for all priceable products
    1b
 };

.portfolio.__calculateSingleTradeGreeksSafely:{[tradeDictionary;marketData;model;config]
    barrierType:.product.getBarrierType tradeDictionary;
    metaFields:`tradeId`underlying`productType`exerciseStyle`optionType`barrierType!(
        tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
        tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType);
    if[not .portfolio.__isGreeksSupported tradeDictionary;
        :metaFields,`delta`gamma`theta`vega`rho`status`statusMessage!(
            0Nf;0Nf;0Nf;0Nf;0Nf;`UNSUPPORTED;
            "Portfolio Greeks v0.8 only supports European vanilla options")];
    functionResult:.[.greeks.calculateGreeks;(tradeDictionary;marketData;model;config);{x}];
    if[10h=type functionResult;
        :metaFields,`delta`gamma`theta`vega`rho`status`statusMessage!(
            0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;functionResult)];
    metaFields,`delta`gamma`theta`vega`rho`status`statusMessage!(
        functionResult[`delta]0;functionResult[`gamma]0;functionResult[`theta]0;
        functionResult[`vega]0;functionResult[`rho]0;`OK;"")
 };

.portfolio.__generateSingleTradeScenarioSafely:{[tradeDictionary;marketData;model;config]
    barrierType:.product.getBarrierType tradeDictionary;
    functionResult:.[.risk.generateScenarioReport;(tradeDictionary;marketData;model;config);{x}];
    if[10h=type functionResult;
        :enlist `tradeId`underlying`productType`exerciseStyle`optionType`barrierType`scenario`spot`riskFreeRate`dividendYield`volatility`unitPrice`unitPnL`notionalPrice`notionalPnL`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            `error;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;functionResult)];
    nRows:count functionResult;
    ([] tradeId:functionResult`tradeId;
        underlying:functionResult`underlying;
        productType:nRows#tradeDictionary`productType;
        exerciseStyle:nRows#tradeDictionary`exerciseStyle;
        optionType:functionResult`optionType;
        barrierType:nRows#barrierType;
        scenario:functionResult`scenario;
        spot:functionResult`spot;
        riskFreeRate:functionResult`riskFreeRate;
        dividendYield:functionResult`dividendYield;
        volatility:functionResult`volatility;
        unitPrice:functionResult`unitPrice;
        unitPnL:functionResult`unitPnL;
        notionalPrice:functionResult`notionalPrice;
        notionalPnL:functionResult`notionalPnL;
        status:nRows#`OK;
        statusMessage:nRows#enlist "")
 };

/ ============================================================================
/ Market data book portfolio functions (v0.10)
/ ============================================================================

.portfolio.priceTradesWithMarketDataBook:{[tradeTable;marketDataBook;model;config]
    .portfolio.__validateTradeTable tradeTable;
    .marketbook.validateMarketDataBook marketDataBook;
    numTrades:count tradeTable;
    resultTable:();
    loopIdx:0;
    while[loopIdx<numTrades;
        currentTrade:tradeTable loopIdx;
        singleResult:.portfolio.__priceSingleTradeFromBookSafely[currentTrade;marketDataBook;model;config];
        resultTable:resultTable,enlist singleResult;
        loopIdx+:1];
    resultTable
 };

.portfolio.generatePortfolioScenarioReportWithMarketDataBook:{[tradeTable;marketDataBook;model;config]
    .portfolio.__validateTradeTable tradeTable;
    .marketbook.validateMarketDataBook marketDataBook;
    numTrades:count tradeTable;
    resultTable:();
    loopIdx:0;
    while[loopIdx<numTrades;
        currentTrade:tradeTable loopIdx;
        singleBlock:.portfolio.__generateScenarioFromBookSafely[currentTrade;marketDataBook;model;config];
        resultTable:resultTable,singleBlock;
        loopIdx+:1];
    resultTable
 };

.portfolio.__priceSingleTradeFromBookSafely:{[tradeDictionary;marketDataBook;model;config]
    barrierType:.product.getBarrierType tradeDictionary;
    functionResult:.[{[td;mdb;mdl;cfg]
        marketDataForTrade:.marketbook.getMarketDataForTrade[mdb;td];
        .engine.priceOption[td;marketDataForTrade;mdl;cfg]
     };(tradeDictionary;marketDataBook;model;config);{x}];
    if[10h=type functionResult;
        :`tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            0Nf;0Nf;config`method;model`modelName;`ERROR;functionResult)];
    `tradeId`underlying`productType`exerciseStyle`optionType`barrierType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
        tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
        tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
        functionResult`unitPrice;functionResult`notionalPrice;
        functionResult`method;functionResult`modelName;`OK;"")
 };

.portfolio.__generateScenarioFromBookSafely:{[tradeDictionary;marketDataBook;model;config]
    barrierType:.product.getBarrierType tradeDictionary;
    functionResult:.[{[td;mdb;mdl;cfg]
        marketDataForTrade:.marketbook.getMarketDataForTrade[mdb;td];
        .risk.generateScenarioReport[td;marketDataForTrade;mdl;cfg]
     };(tradeDictionary;marketDataBook;model;config);{x}];
    if[10h=type functionResult;
        :enlist `tradeId`underlying`productType`exerciseStyle`optionType`barrierType`scenario`spot`riskFreeRate`dividendYield`volatility`unitPrice`unitPnL`notionalPrice`notionalPnL`status`statusMessage!(
            tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
            tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType;
            `error;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;functionResult)];
    nRows:count functionResult;
    ([] tradeId:functionResult`tradeId;
        underlying:functionResult`underlying;
        productType:nRows#tradeDictionary`productType;
        exerciseStyle:nRows#tradeDictionary`exerciseStyle;
        optionType:functionResult`optionType;
        barrierType:nRows#barrierType;
        scenario:functionResult`scenario;
        spot:functionResult`spot;
        riskFreeRate:functionResult`riskFreeRate;
        dividendYield:functionResult`dividendYield;
        volatility:functionResult`volatility;
        unitPrice:functionResult`unitPrice;
        unitPnL:functionResult`unitPnL;
        notionalPrice:functionResult`notionalPrice;
        notionalPnL:functionResult`notionalPnL;
        status:nRows#`OK;
        statusMessage:nRows#enlist "")
 };
