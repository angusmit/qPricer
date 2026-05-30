\l core/init.q

-1 "=============================================================================";
-1 " qFDM v0.30 Daily Risk Run";
-1 "=============================================================================";
-1 "";

/ Build sample data
spotTable:([] underlying:`AAPL`MSFT; spot:100 250f);
volTable:([] underlying:`AAPL`MSFT; volatility:0.2 0.25);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:`AAPL`MSFT; dividendYield:0 0.01f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
prevMktBook:mktBook;
shockTable:.histscen.syntheticShockTable[`AAPL`MSFT];

tradeTable:();
tradeTable:tradeTable,enlist `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`bookName`barrierType`barrierLevel`rebate!(1;`AAPL;`equityOption;`european;`call;100f;1f;100000f;`equities;`none;0Nf;0f);
tradeTable:tradeTable,enlist `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`bookName`barrierType`barrierLevel`rebate!(2;`MSFT;`equityOption;`european;`call;250f;1f;50000f;`equities;`none;0Nf;0f);

limitTable:();
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`VaR95;500000f;0.8;1.0;`lessThan;1b);
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(2;`portfolio;`ALL;`ES95;500000f;0.8;1.0;`lessThan;1b);

configDict:`confidenceLevels`mcConfig!(0.95 0.99;`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;25;42;0b;0b;0.95));

/ Run
dailyResult:.dailyrisk.runDailyRisk[tradeTable;mktBook;prevMktBook;shockTable;limitTable;configDict];

/ Dashboard
dashSummary:dailyResult`dashboardSummary;
-1 "--- Dashboard Summary ---";
-1 "Trades: ",string[dashSummary`tradeCount],"  PricedOK: ",string[dashSummary`pricedOkCount],"  TotalPV: ",string dashSummary`totalPV;
-1 "VaR95: ",string[dashSummary`var95],"  ES95: ",string dashSummary`es95;
-1 "HistVaR95: ",string[dashSummary`historicalVar95],"  WorstHistLoss: ",string dashSummary`worstHistoricalLoss;
-1 "Limits OK: ",string[dashSummary`limitOkCount],"  Warnings: ",string[dashSummary`limitWarningCount],"  Breaches: ",string dashSummary`limitBreachCount;
-1 "ModelChecks Passed: ",string[dashSummary`modelCheckPassedCount],"  Failed: ",string dashSummary`modelCheckFailedCount;
-1 "Overall: ",string dashSummary`overallStatus;
-1 "";

/ Run status
-1 "--- Run Status ---";
runStatusTable:dailyResult`runStatus;
rsIdx:0;
while[rsIdx<count runStatusTable;
    rsRow:runStatusTable rsIdx;
    -1 "  ",string[rsRow`blockName]," ",string[rsRow`status]," rows=",string[rsRow`rowCount]," err=",string[rsRow`errorCount]," ",string[rsRow`elapsedMs],"ms";
    rsIdx+:1];
-1 "";

/ Limit breaches
breaches:.limits.limitBreaches dailyResult`limitCheckResult;
if[0<count breaches;
    -1 "--- Limit Breaches ---";
    bIdx:0;
    while[bIdx<count breaches;
        bRow:breaches bIdx;
        -1 "  ",string[bRow`metricName]," ",string[bRow`metricValue]," vs limit ",string[bRow`limitValue]," util=",string bRow`utilisation;
        bIdx+:1];
    -1 ""];

/ Worst events
worstEvents:dailyResult`worstHistoricalEvents;
if[0<count worstEvents;
    -1 "--- Worst Historical Events ---";
    wIdx:0;
    wCnt:3&count worstEvents;
    while[wIdx<wCnt;
        wRow:worstEvents wIdx;
        -1 "  #",string[wRow`rank]," ",string[wRow`eventName]," pnl=",string[wRow`pnl]," loss=",string wRow`loss;
        wIdx+:1];
    -1 ""];

-1 "=============================================================================";
-1 " Daily risk run completed: ",string dashSummary`overallStatus;
-1 "=============================================================================";
