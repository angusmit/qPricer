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
pnlDist:.replay.historicalPnlDistribution replayResult;
scenCount:count .histscen.scenarioKeys shockTable;
.testutil.assertTrue[scenCount=count pnlDist;"one row per scenario"];
pnlVec:pnlDist`pnl;
lossVec:pnlDist`loss;
.testutil.assertTrue[all lossVec=(neg pnlVec);"loss = -pnl"];
/ Crash scenarios should produce net losses (negative PnL) for call portfolio
.testutil.assertTrue[any pnlVec<0f;"at least one losing scenario"];
.testutil.assertTrue[any lossVec>0f;"at least one positive loss"];
-1 "PASS test_historical_replay_distribution: scenarios=",string count pnlDist;
