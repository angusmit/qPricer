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

.portfolio.__isGreeksSupported:{[tradeDictionary]
    if[not tradeDictionary[`exerciseStyle]~`european; :0b];
    if[.product.isBarrierOption tradeDictionary; :0b];
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
