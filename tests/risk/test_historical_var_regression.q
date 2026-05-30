\l core/init.q
/ Deterministic sign convention regression test
/ Synthetic replay result with known PnLs
syntheticReplay:();
syntheticReplay:syntheticReplay,enlist `scenarioDate`eventName`tradeId`bookName`underlying`productType`basePV`shockedPV`historicalPnl`status`errorMessage!(2020.03.16;`covidCrash;1;`equities;`AAPL;`equityOption;1000f;900f;-100f;`OK;"");
syntheticReplay:syntheticReplay,enlist `scenarioDate`eventName`tradeId`bookName`underlying`productType`basePV`shockedPV`historicalPnl`status`errorMessage!(2021.01.04;`rally;1;`equities;`AAPL;`equityOption;1000f;1050f;50f;`OK;"");
syntheticReplay:syntheticReplay,enlist `scenarioDate`eventName`tradeId`bookName`underlying`productType`basePV`shockedPV`historicalPnl`status`errorMessage!(2022.06.13;`rateShock;1;`equities;`AAPL;`equityOption;1000f;980f;-20f;`OK;"");

pnlDist:.replay.historicalPnlDistribution syntheticReplay;
pnlVec:pnlDist`pnl;
lossVec:pnlDist`loss;
/ Sign checks
.testutil.assertTrue[(pnlDist 0)[`loss]=100f;"covidCrash loss=100"];
.testutil.assertTrue[(pnlDist 1)[`loss]=-50f;"rally loss=-50"];
.testutil.assertTrue[(pnlDist 2)[`loss]=20f;"rateShock loss=20"];

/ VaR should be positive (portfolio has losses)
varRep:.replay.historicalVarReport[syntheticReplay;enlist 0.5];
varRow:varRep 0;
.testutil.assertTrue[varRow[`valueAtRisk]>0f;"VaR positive"];
.testutil.assertTrue[varRow[`expectedShortfall]>0f;"ES positive"];

/ Worst event
worstRows:.replay.worstHistoricalEvents[syntheticReplay;3];
.testutil.assertTrue[(worstRows 0)[`loss]=100f;"worst loss = 100"];
.testutil.assertTrue[(worstRows 0)[`eventName]=`covidCrash;"worst event is covidCrash"];

-1 "PASS test_historical_var_regression: VaR=",string[varRow`valueAtRisk],", worstLoss=",string (worstRows 0)`loss;
