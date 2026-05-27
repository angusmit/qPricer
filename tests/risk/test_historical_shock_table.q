\l lib/init.q
shockTable:.histscen.syntheticShockTable[`AAPL`MSFT];
.testutil.assertTrue[(count shockTable)>0;"synthetic table non-empty"];
.histscen.validateHistoricalShockTable shockTable;
scenKeys:.histscen.scenarioKeys shockTable;
.testutil.assertTrue[(count scenKeys)>0;"scenario keys non-empty"];
shockRows:.histscen.shocksForScenario[shockTable;2020.03.16;`covidCrash];
.testutil.assertTrue[(count shockRows)>0;"shock rows for COVID"];
spotTable:([] underlying:`AAPL`MSFT; spot:100 250f);
volTable:([] underlying:`AAPL`MSFT; volatility:0.2 0.25);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:`AAPL`MSFT; dividendYield:0 0.01f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
shockedBook:.histscen.applyHistoricalShock[mktBook;shockRows;()!()];
shockedSpot:.marketbook.getSpot[shockedBook;`AAPL];
.testutil.assertTrue[shockedSpot<100f;"AAPL spot shocked down"];
-1 "PASS test_historical_shock_table: shockedSpot=",string shockedSpot;
