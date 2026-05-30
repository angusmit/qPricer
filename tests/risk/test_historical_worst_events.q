\l core/init.q
spotTable:([] underlying:`AAPL`MSFT; spot:100 250f);
volTable:([] underlying:`AAPL`MSFT; volatility:0.2 0.25);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:`AAPL`MSFT; dividendYield:0 0.01f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
shockTable:.histscen.syntheticShockTable[`AAPL`MSFT];
tradeTable:();
tradeTable:tradeTable,enlist `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`bookName`barrierType`barrierLevel`rebate!(1;`AAPL;`equityOption;`european;`call;100f;1f;100000f;`equities;`none;0Nf;0f);
tradeTable:tradeTable,enlist `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`bookName`barrierType`barrierLevel`rebate!(2;`MSFT;`equityOption;`european;`call;250f;1f;50000f;`equities;`none;0Nf;0f);
replayResult:.replay.replayHistoricalScenarios[tradeTable;mktBook;shockTable;()!()];
worstRows:.replay.worstHistoricalEvents[replayResult;2];
.testutil.assertTrue[2=count worstRows;"top 2 events"];
.testutil.assertTrue[(worstRows 0)[`rank]=1;"rank 1"];
.testutil.assertTrue[(worstRows 0)[`loss]>=(worstRows 1)`loss;"sorted by loss"];
.testutil.assertTrue[(worstRows 0)[`loss]>0f;"worst loss is positive"];
.testutil.assertTrue[`scenarioDate in key worstRows 0;"has scenarioDate"];
-1 "PASS test_historical_worst_events: worstLoss=",string (worstRows 0)`loss;
