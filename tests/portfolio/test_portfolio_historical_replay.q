\l lib/init.q
spotTable:([] underlying:`AAPL`MSFT; spot:100 250f);
volTable:([] underlying:`AAPL`MSFT; volatility:0.2 0.25);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:`AAPL`MSFT; dividendYield:0 0.01f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
shockTable:.histscen.syntheticShockTable[`AAPL`MSFT];
tradeTable:();
tradeTable:tradeTable,enlist `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`bookName`barrierType`barrierLevel`rebate!(1;`AAPL;`equityOption;`european;`call;100f;1f;100000f;`equities;`none;0Nf;0f);
tradeTable:tradeTable,enlist `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`bookName`barrierType`barrierLevel`rebate!(2;`MSFT;`equityOption;`european;`put;250f;1f;50000f;`equities;`none;0Nf;0f);
replayResult:.replay.replayHistoricalScenarios[tradeTable;mktBook;shockTable;()!()];
summaryResult:.replay.replaySummary replayResult;
.testutil.assertTrue[summaryResult[`scenarioCount]>0;"scenarios ran"];
.testutil.assertTrue[summaryResult[`okRows]>0;"some OK"];
varRep:.replay.historicalVarReport[replayResult;0.5 0.9];
.testutil.assertTrue[2=count varRep;"VaR report rows"];
worstEvents:.replay.worstHistoricalEvents[replayResult;3];
.testutil.assertTrue[(count worstEvents)>0;"worst events exist"];
-1 "PASS test_portfolio_historical_replay: scenarios=",string[summaryResult`scenarioCount],", okRows=",string summaryResult`okRows;
