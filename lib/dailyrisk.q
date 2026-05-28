/ dailyrisk.q - end-to-end daily risk workflow (v0.30)

.dailyrisk.__blockResult:{[blockName;statusSym;rowCnt;errCnt;errorMsg]
    `blockName`status`rowCount`errorCount`warningCount`elapsedMs`errorMessage!(
        blockName;statusSym;rowCnt;errCnt;0;0;"",errorMsg)
 };

.dailyrisk.validateDailyRiskInputs:{[tradeTable;marketDataBook;previousMarketDataBook;historicalShockTable;limitTable;configDict]
    if[0=count tradeTable; '"Trade table is empty"];
    if[not 99h=type marketDataBook; '"Invalid marketDataBook"];
 };

.dailyrisk.runPricingBlock:{[tradeTable;marketDataBook;configDict]
    bsModel:.model.createBlackScholesModel[];
    fdmConfig:.replay.__stableConfig[];
    priceFn:{.portfolio.priceTradesWithMarketDataBook[x 0;x 1;x 2;x 3]};
    pricingResult:@[priceFn;(tradeTable;marketDataBook;bsModel;fdmConfig);{x}];
    if[10h=type pricingResult;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`pricing;`ERROR;0;0;pricingResult]))];
    okCnt:sum (pricingResult`status)=`OK;
    errCnt:(count pricingResult)-okCnt;
    `result`runStatus!(pricingResult;.dailyrisk.__blockResult[`pricing;`OK;count pricingResult;errCnt;""])
 };

.dailyrisk.runGreekBlock:{[tradeTable;marketDataBook;configDict]
    bsModel:.model.createBlackScholesModel[];
    fdmConfig:.replay.__stableConfig[];
    resultRows:();
    tIdx:0;
    while[tIdx<count tradeTable;
          tradeRow:tradeTable tIdx;
          greekFn:{
              mkt:.marketbook.getMarketDataForTrade[x 1;x 0];
              .greeks.calculateGreeks[x 0;mkt;x 2;x 3]};
          greekResult:@[greekFn;(tradeRow;marketDataBook;bsModel;fdmConfig);{x}];
          if[not 10h=type greekResult;
             greekRow:greekResult 0;
             resultRows:resultRows,enlist greekRow];
          tIdx+:1];
    `result`runStatus!(resultRows;.dailyrisk.__blockResult[`greeks;`OK;count resultRows;0;""])
 };

.dailyrisk.runScenarioBlock:{[tradeTable;marketDataBook;configDict]
    bsModel:.model.createBlackScholesModel[];
    fdmConfig:.replay.__stableConfig[];
    scenFn:{.portfolio.generatePortfolioScenarioReportWithMarketDataBook[x 0;x 1;x 2;x 3]};
    scenResult:@[scenFn;(tradeTable;marketDataBook;bsModel;fdmConfig);{x}];
    if[10h=type scenResult;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`scenario;`ERROR;0;0;scenResult]))];
    `result`runStatus!(scenResult;.dailyrisk.__blockResult[`scenario;`OK;count scenResult;0;""])
 };

.dailyrisk.runPnlExplainBlock:{[tradeTable;marketDataBook;previousMarketDataBook;configDict]
    if[(::)~previousMarketDataBook;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`pnlExplain;`skipped;0;0;"No previous market data"]))];
    pnlFn:{.pnl.explainPortfolio[x 0;x 1;x 2;x 3]};
    pnlResult:@[pnlFn;(tradeTable;previousMarketDataBook;marketDataBook;configDict);{x}];
    if[10h=type pnlResult;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`pnlExplain;`ERROR;0;0;pnlResult]))];
    `result`runStatus!(pnlResult;.dailyrisk.__blockResult[`pnlExplain;`OK;count pnlResult;0;""])
 };

.dailyrisk.runVarBlock:{[scenarioResult;confidenceLevels]
    if[0=count scenarioResult;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`var;`ERROR;0;0;"No scenario data"]))];
    / Map scenario report columns to VaR engine format if needed
    firstRow:scenarioResult 0;
    rowKeys:key firstRow;
    mappedResult:scenarioResult;
    if[(`scenario in rowKeys) and not `scenarioName in rowKeys;
       mapped:();
       mIdx:0;
       while[mIdx<count scenarioResult;
             srcRow:scenarioResult mIdx;
             pnlVal:$[`scenarioPnl in rowKeys;srcRow`scenarioPnl;
                      $[`notionalPnL in rowKeys;srcRow`notionalPnL;0Nf]];
             mapped:mapped,enlist `scenarioName`scenarioPnl!(srcRow`scenario;pnlVal);
             mIdx+:1];
       mappedResult:mapped];
    varFn:{.var.varReportFromScenarioResult[x 0;x 1]};
    varReport:@[varFn;(mappedResult;confidenceLevels);{x}];
    if[10h=type varReport;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`var;`ERROR;0;0;varReport]))];
    `result`runStatus!(varReport;.dailyrisk.__blockResult[`var;`OK;count varReport;0;""])
 };

.dailyrisk.runHistoricalReplayBlock:{[tradeTable;marketDataBook;historicalShockTable;configDict]
    if[(::)~historicalShockTable;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`historicalReplay;`skipped;0;0;"No historical shock table"]))];
    replayFn:{.replay.replayHistoricalScenarios[x 0;x 1;x 2;x 3]};
    replayResult:@[replayFn;(tradeTable;marketDataBook;historicalShockTable;configDict);{x}];
    if[10h=type replayResult;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`historicalReplay;`ERROR;0;0;replayResult]))];
    `result`runStatus!(replayResult;.dailyrisk.__blockResult[`historicalReplay;`OK;count replayResult;0;""])
 };

.dailyrisk.runLimitBlock:{[varReport;greekReport;pnlReport;historicalReplayResult;limitTable;configDict]
    if[(::)~limitTable;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`limits;`skipped;0;0;"No limit table"]))];
    if[0=count limitTable;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`limits;`skipped;0;0;"Empty limit table"]))];
    histReport:$[(::)~historicalReplayResult;(::);
                 @[{.replay.worstHistoricalEvents[x 0;x 1]};(historicalReplayResult;5);{(::)}]];
    metricTable:.limits.buildMetricTable[varReport;greekReport;pnlReport;histReport];
    if[0=count metricTable;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`limits;`ERROR;0;0;"Empty metric table"]))];
    limitFn:{.limits.evaluateLimits[x 0;x 1]};
    limitResult:@[limitFn;(limitTable;metricTable);{x}];
    if[10h=type limitResult;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`limits;`ERROR;0;0;limitResult]))];
    `result`runStatus!(limitResult;.dailyrisk.__blockResult[`limits;`OK;count limitResult;0;""])
 };

.dailyrisk.runModelCheckBlock:{[configDict]
    mcFn:{.modelcheck.runCoreModelLimitChecks x};
    checkResult:@[mcFn;configDict;{x}];
    if[10h=type checkResult;
       :(`result`runStatus!(();.dailyrisk.__blockResult[`modelCheck;`ERROR;0;0;checkResult]))];
    `result`runStatus!(checkResult;.dailyrisk.__blockResult[`modelCheck;`OK;count checkResult;0;""])
 };

.dailyrisk.runDailyRisk:{[tradeTable;marketDataBook;previousMarketDataBook;historicalShockTable;limitTable;configDict]
    .dailyrisk.validateDailyRiskInputs[tradeTable;marketDataBook;previousMarketDataBook;historicalShockTable;limitTable;configDict];
    confidenceLevels:$[0<count configDict;
                       $[`confidenceLevels in key configDict;configDict`confidenceLevels;0.95 0.99];
                       0.95 0.99];
    statusRows:();
    / 1. Pricing
    pricingBlock:.dailyrisk.runPricingBlock[tradeTable;marketDataBook;configDict];
    pricingResult:pricingBlock`result;
    statusRows:statusRows,enlist pricingBlock`runStatus;
    / 2. Greeks
    greekBlock:.dailyrisk.runGreekBlock[tradeTable;marketDataBook;configDict];
    greekResult:greekBlock`result;
    statusRows:statusRows,enlist greekBlock`runStatus;
    / 3. Scenarios
    scenarioBlock:.dailyrisk.runScenarioBlock[tradeTable;marketDataBook;configDict];
    scenarioResult:scenarioBlock`result;
    statusRows:statusRows,enlist scenarioBlock`runStatus;
    / 4. PnL explain
    pnlBlock:.dailyrisk.runPnlExplainBlock[tradeTable;marketDataBook;previousMarketDataBook;configDict];
    pnlExplainResult:pnlBlock`result;
    statusRows:statusRows,enlist pnlBlock`runStatus;
    / 5. VaR
    varBlock:.dailyrisk.runVarBlock[scenarioResult;confidenceLevels];
    varReport:varBlock`result;
    statusRows:statusRows,enlist varBlock`runStatus;
    / 6. Historical replay
    histBlock:.dailyrisk.runHistoricalReplayBlock[tradeTable;marketDataBook;historicalShockTable;configDict];
    histReplayResult:histBlock`result;
    statusRows:statusRows,enlist histBlock`runStatus;
    / Historical PnL distribution / VaR / worst events
    histPnlDist:$[0<count histReplayResult;
                  @[.replay.historicalPnlDistribution;histReplayResult;{()}];()];
    histVarReport:$[0<count histReplayResult;
                    @[{.replay.historicalVarReport[x 0;x 1]};(histReplayResult;confidenceLevels);{()}];()];
    worstEvents:$[0<count histReplayResult;
                  @[{.replay.worstHistoricalEvents[x 0;x 1]};(histReplayResult;5);{()}];()];
    statusRows:statusRows,enlist .dailyrisk.__blockResult[`historicalVar;$[0<count histVarReport;`OK;`skipped];count histVarReport;0;""];
    / 7. Limits
    limitBlock:.dailyrisk.runLimitBlock[varReport;greekResult;();histReplayResult;limitTable;configDict];
    limitCheckResult:limitBlock`result;
    statusRows:statusRows,enlist limitBlock`runStatus;
    limitDash:$[0<count limitCheckResult;.limitreport.limitDashboard limitCheckResult;
                `totalLimits`okCount`warningCount`breachCount`errorCount`worstSeverity`maxUtilisation`status`errorMessage!(0;0;0;0;0;`OK;0Nf;`OK;"")];
    / 8. Model checks
    mcBlock:.dailyrisk.runModelCheckBlock[configDict];
    modelCheckResult:mcBlock`result;
    statusRows:statusRows,enlist mcBlock`runStatus;
    mcSummary:$[0<count modelCheckResult;.limitcheck.summary modelCheckResult;
                `checkCount`passedCount`failedCount`maxAbsoluteDifference`maxRelativeDifference`status!(0;0;0;0Nf;0Nf;`OK)];
    / 9. Dashboard
    dailyResult:`pricingResult`greekResult`scenarioResult`pnlExplainResult`varReport`historicalReplayResult`historicalPnlDistribution`historicalVarReport`worstHistoricalEvents`limitCheckResult`limitDashboard`modelCheckResult`modelCheckSummary`runStatus!(
        pricingResult;greekResult;scenarioResult;pnlExplainResult;varReport;histReplayResult;histPnlDist;histVarReport;worstEvents;limitCheckResult;limitDash;modelCheckResult;mcSummary;statusRows);
    dashSummary:.dashboard.dashboardSummary dailyResult;
    dailyResult[`dashboardSummary]:dashSummary;
    statusRows:statusRows,enlist .dailyrisk.__blockResult[`dashboard;`OK;1;0;""];
    dailyResult[`runStatus]:statusRows;
    dailyResult
 };
