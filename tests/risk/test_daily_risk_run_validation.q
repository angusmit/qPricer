\l core/init.q
/ Empty trade table
emptyResult:@[{.dailyrisk.runDailyRisk[x 0;x 1;x 2;x 3;x 4;x 5]};(();`a`b!(1;2);(::);(::);(::);()!());{`ERROR}];
.testutil.assertTrue[emptyResult~`ERROR;"empty trades fails"];

/ Missing previous market data -> PnL explain skipped
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
noHistResult:.dailyrisk.runDailyRisk[tradeTable;mktBook;(::);(::);(::);configDict];
runStatusTable:noHistResult`runStatus;
pnlStatus:runStatusTable where (runStatusTable`blockName)=`pnlExplain;
.testutil.assertTrue[1=count pnlStatus;"pnl explain row exists"];
.testutil.assertTrue[(pnlStatus 0)[`status]=`skipped;"pnl explain skipped"];

histStatus:runStatusTable where (runStatusTable`blockName)=`historicalReplay;
.testutil.assertTrue[(histStatus 0)[`status]=`skipped;"historical replay skipped"];

limitStatus:runStatusTable where (runStatusTable`blockName)=`limits;
.testutil.assertTrue[(limitStatus 0)[`status]=`skipped;"limits skipped"];

-1 "PASS test_daily_risk_run_validation";
