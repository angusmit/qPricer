/ perfdiag.q - performance diagnostics: timing, scaling, bottleneck detection

/ --- Public: timing ---

.perfdiag.timeStep:{[stepLabel;functionHandle;argumentList]
    startTime:.z.p;
    functionResult:.[functionHandle;argumentList;{x}];
    endTime:.z.p;
    elapsedMs:(`long$endTime-startTime)%1000000;
    statusSym:`OK;
    errorMsg:"";
    if[10h=type functionResult; statusSym:`ERROR; errorMsg:functionResult];
    `stepName`elapsedMs`status`errorMessage!(stepLabel;elapsedMs;statusSym;errorMsg)
 };

.perfdiag.timeBatchRun:{[tradeTable;marketDataBook;previousMarketDataBook;configDict]
    pricingModel:configDict`model;
    fdmConfig:configDict`fdmConfig;
    timeStepYears:configDict`timeStepYears;
    bookLabel:configDict`bookName;
    timingRows:();
    pricingTimed:.perfdiag.timeStep["pricing";
        .portfolio.priceTradesWithMarketDataBook;(tradeTable;marketDataBook;pricingModel;fdmConfig)];
    timingRows:timingRows,enlist pricingTimed;
    greeksTimed:.perfdiag.timeStep["greeks";
        .batch.__calculateGreeksFromBook;(tradeTable;marketDataBook;pricingModel;fdmConfig)];
    timingRows:timingRows,enlist greeksTimed;
    scenarioTimed:.perfdiag.timeStep["scenarios";
        .portfolio.generatePortfolioScenarioReportWithMarketDataBook;(tradeTable;marketDataBook;pricingModel;fdmConfig)];
    timingRows:timingRows,enlist scenarioTimed;
    pnlConfig:`model`fdmConfig`timeStepYears`bookName!(pricingModel;fdmConfig;timeStepYears;bookLabel);
    pnlTimed:.perfdiag.timeStep["pnlExplain";
        .pnl.explainPortfolio;(tradeTable;previousMarketDataBook;marketDataBook;pnlConfig)];
    timingRows:timingRows,enlist pnlTimed;
    pricingResult:.portfolio.priceTradesWithMarketDataBook[tradeTable;marketDataBook;pricingModel;fdmConfig];
    reportTimed:.perfdiag.timeStep["reports";
        .report.portfolioSummary;enlist pricingResult];
    timingRows:timingRows,enlist reportTimed;
    timingRows
 };

/ --- Public: scaling tests ---

.perfdiag.portfolioScalingTest:{[tradeCounts;symbolCountValue;configDict;modeSym]
    resultList:();
    loopIdx:0;
    while[loopIdx<count tradeCounts;
        tradeCountValue:tradeCounts loopIdx;
        stressResult:.stress.runPortfolioStress[tradeCountValue;symbolCountValue;configDict;modeSym];
        scalingRow:`tradeCount`symbolCount`elapsedMs`okRows`errorRows!(
            tradeCountValue;symbolCountValue;stressResult`elapsedMs;stressResult`okRows;stressResult`errorRows);
        resultList:resultList,enlist scalingRow;
        loopIdx+:1];
    resultList
 };

.perfdiag.scenarioScalingTest:{[symbolCountValue;tradeCounts;configDict]
    pricingModel:configDict`model;
    fdmConfig:configDict`fdmConfig;
    scenariosPerTrade:9;
    resultList:();
    loopIdx:0;
    while[loopIdx<count tradeCounts;
        currentTradeCount:tradeCounts loopIdx;
        marketDataBook:.stress.generateMarketDataBook[symbolCountValue;2025.01.01];
        symbolList:marketDataBook[`spotTable]`underlying;
        tradeTable:.stress.generateSupportedTradeTable[currentTradeCount;symbolList;2025.01.01];
        expectedRowCount:currentTradeCount*scenariosPerTrade;
        startTime:.z.p;
        scenarioResult:.portfolio.generatePortfolioScenarioReportWithMarketDataBook[tradeTable;marketDataBook;pricingModel;fdmConfig];
        elapsedMs:(`long$.z.p-startTime)%1000000;
        actualRowCount:count scenarioResult;
        droppedRowCount:expectedRowCount-actualRowCount;
        okRowCount:sum scenarioResult[`status]=`OK;
        errorRowCount:actualRowCount-okRowCount;
        scalingRow:`tradeCount`scenarioCount`expectedRows`actualRows`droppedRows`elapsedMs`okRows`errorRows!(
            currentTradeCount;scenariosPerTrade;expectedRowCount;actualRowCount;droppedRowCount;
            elapsedMs;okRowCount;errorRowCount);
        resultList:resultList,enlist scalingRow;
        loopIdx+:1];
    resultList
 };

.perfdiag.gridScalingTest:{[tradeRow;marketDataBook;pricingModel;configList]
    / Get market data and BS benchmark once
    marketDataForTrade:.marketbook.getMarketDataForTrade[marketDataBook;tradeRow];
    benchmarkPrice:.validation.blackScholesClosedForm[tradeRow`optionType;
        marketDataForTrade`spot;tradeRow`strike;tradeRow`expiry;
        marketDataForTrade`riskFreeRate;marketDataForTrade`dividendYield;
        marketDataForTrade`volatility];
    resultList:();
    loopIdx:0;
    while[loopIdx<count configList;
        currentConfig:configList loopIdx;
        startTime:.z.p;
        priceResult:.[.engine.priceOption;(tradeRow;marketDataForTrade;pricingModel;currentConfig);{x}];
        elapsedMs:(`long$.z.p-startTime)%1000000;
        fdmPrice:0Nf; absError:0Nf; statusSym:`ERROR; errorMsg:"";
        if[99h=type priceResult;
            fdmPrice:priceResult`unitPrice;
            absError:abs fdmPrice-benchmarkPrice;
            statusSym:`OK];
        if[10h=type priceResult; errorMsg:priceResult];
        gridRow:`numberOfSpotSteps`numberOfTimeSteps`elapsedMs`unitPrice`benchmarkPrice`absoluteError`status`errorMessage!(
            currentConfig`numberOfSpotSteps;currentConfig`numberOfTimeSteps;elapsedMs;
            fdmPrice;benchmarkPrice;absError;statusSym;errorMsg);
        resultList:resultList,enlist gridRow;
        loopIdx+:1];
    resultList
 };

/ --- Public: analysis ---

.perfdiag.summariseTimings:{[timingTable]
    totalMs:sum timingTable`elapsedMs;
    stepCount:count timingTable;
    okSteps:sum timingTable[`status]=`OK;
    errorSteps:stepCount-okSteps;
    `totalMs`stepCount`okSteps`errorSteps!(totalMs;stepCount;okSteps;errorSteps)
 };

.perfdiag.detectSlowSteps:{[timingTable;thresholdMs]
    timingTable where timingTable[`elapsedMs]>thresholdMs
 };
