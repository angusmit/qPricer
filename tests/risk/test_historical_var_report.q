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
varRep:.replay.historicalVarReport[replayResult;0.5 0.9];
.testutil.assertTrue[2=count varRep;"2 confidence levels"];
firstRow:varRep 0;
.testutil.assertTrue[not null firstRow`valueAtRisk;"VaR finite"];
.testutil.assertTrue[not null firstRow`expectedShortfall;"ES finite"];
.testutil.assertTrue[firstRow[`status]~`OK;"status OK"];
-1 "PASS test_historical_var_report: VaR=",string firstRow`valueAtRisk;
