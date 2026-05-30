\l core/init.q
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
dailyResult:.dailyrisk.runDailyRisk[tradeTable;mktBook;prevMktBook;shockTable;limitTable;configDict];
.testutil.assertTrue[`pricingResult in key dailyResult;"has pricingResult"];
.testutil.assertTrue[`varReport in key dailyResult;"has varReport"];
.testutil.assertTrue[`historicalReplayResult in key dailyResult;"has historicalReplayResult"];
.testutil.assertTrue[`limitCheckResult in key dailyResult;"has limitCheckResult"];
.testutil.assertTrue[`dashboardSummary in key dailyResult;"has dashboardSummary"];
.testutil.assertTrue[`runStatus in key dailyResult;"has runStatus"];

pricingRows:dailyResult`pricingResult;
.testutil.assertTrue[(count pricingRows)>0;"pricing has rows"];

runStatusTable:dailyResult`runStatus;
.testutil.assertTrue[(count runStatusTable)>0;"runStatus has rows"];

dashSummary:dailyResult`dashboardSummary;
.testutil.assertTrue[dashSummary[`tradeCount]=2;"2 trades"];

-1 "PASS test_daily_risk_run: trades=",string[dashSummary`tradeCount],", overallStatus=",string dashSummary`overallStatus;
