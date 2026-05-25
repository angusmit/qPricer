/ batch.q - daily pricing/risk batch orchestration

/ --- Public ---

.batch.runDailyPricing:{[tradeTable;marketDataBook;previousMarketDataBook;configDict]
    pricingModel:configDict`model;
    fdmConfig:configDict`fdmConfig;
    timeStepYears:configDict`timeStepYears;
    bookLabel:configDict`bookName;
    valuationDate:configDict`valuationDate;
    runLabel:configDict`runLabel;
    / Price
    pricingResult:.portfolio.priceTradesWithMarketDataBook[tradeTable;marketDataBook;pricingModel;fdmConfig];
    / Greeks (per-trade from book)
    greekResult:.batch.__calculateGreeksFromBook[tradeTable;marketDataBook;pricingModel;fdmConfig];
    / Scenarios
    scenarioResult:.portfolio.generatePortfolioScenarioReportWithMarketDataBook[tradeTable;marketDataBook;pricingModel;fdmConfig];
    / PnL explain
    pnlConfig:`model`fdmConfig`timeStepYears`bookName!(pricingModel;fdmConfig;timeStepYears;bookLabel);
    pnlExplainResult:.pnl.explainPortfolio[tradeTable;previousMarketDataBook;marketDataBook;pnlConfig];
    / Summaries
    portfolioSummary:.report.portfolioSummary pricingResult;
    riskSummary:.report.riskSummary greekResult;
    scenarioSummary:.report.scenarioSummary scenarioResult;
    / Build partial result for audit
    partialResult:`pricingResult`greekResult`scenarioResult`pnlExplainResult`portfolioSummary`riskSummary`scenarioSummary!(
        pricingResult;greekResult;scenarioResult;pnlExplainResult;portfolioSummary;riskSummary;scenarioSummary);
    errorSummary:.audit.errorAudit[pricingResult;scenarioResult;pnlExplainResult];
    auditRecord:.audit.createAuditRecord[partialResult;valuationDate;runLabel];
    partialResult,`errorSummary`auditRecord!(errorSummary;auditRecord)
 };

.batch.writeDailyReports:{[runResult;outputDirectory;valuationDate]
    dateStr:valuationDate;
    if[not 10h=type valuationDate; dateStr:string valuationDate];
    .report.exportCsv[runResult`pricingResult;outputDirectory,"/pricing_",dateStr,".csv"];
    .report.exportCsv[runResult`scenarioResult;outputDirectory,"/scenarios_",dateStr,".csv"];
    .report.exportCsv[runResult`pnlExplainResult;outputDirectory,"/pnl_explain_",dateStr,".csv"];
    -1 "Reports written to ",outputDirectory;
 };

/ --- Internal: Greeks with market data book ---

.batch.__calculateGreeksFromBook:{[tradeTable;marketDataBook;pricingModel;fdmConfig]
    numTrades:count tradeTable;
    resultList:();
    loopIdx:0;
    while[loopIdx<numTrades;
          currentTrade:tradeTable loopIdx;
          singleResult:.batch.__greeksSingleFromBook[currentTrade;marketDataBook;pricingModel;fdmConfig];
          resultList:resultList,enlist singleResult;
          loopIdx+:1];
    resultList
 };

.batch.__greeksSingleFromBook:{[tradeDictionary;marketDataBook;pricingModel;fdmConfig]
    barrierType:.product.getBarrierType tradeDictionary;
    metaFields:`tradeId`underlying`productType`exerciseStyle`optionType`barrierType!(
        tradeDictionary`tradeId;tradeDictionary`underlying;tradeDictionary`productType;
        tradeDictionary`exerciseStyle;tradeDictionary`optionType;barrierType);
    unsupportedRow:metaFields,`delta`gamma`theta`vega`rho`status`statusMessage!(
        0Nf;0Nf;0Nf;0Nf;0Nf;`UNSUPPORTED;"Greeks only supported for European vanilla");
    if[not tradeDictionary[`exerciseStyle]~`european; :unsupportedRow];
    if[.product.isBarrierOption tradeDictionary; :unsupportedRow];
    functionResult:.[{[td;mdb;mdl;cfg]
                         mktData:.marketbook.getMarketDataForTrade[mdb;td];
                         .greeks.calculateGreeks[td;mktData;mdl;cfg]
                     };(tradeDictionary;marketDataBook;pricingModel;fdmConfig);{x}];
    if[10h=type functionResult;
       :metaFields,`delta`gamma`theta`vega`rho`status`statusMessage!(0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;functionResult)];
    metaFields,`delta`gamma`theta`vega`rho`status`statusMessage!(
        functionResult[`delta]0;functionResult[`gamma]0;functionResult[`theta]0;
        functionResult[`vega]0;functionResult[`rho]0;`OK;"")
 };
