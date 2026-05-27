\l lib/init.q
/ Synthetic scenario result: 2 trades x 3 scenarios
scenarioResult:();
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`spotUp;5f);
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`spotUp;3f);
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`spotDown;-8f);
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`spotDown;-2f);
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`base;1f);
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`base;0.5);

pnlDist:.riskdist.buildPnlDistribution scenarioResult;
.testutil.assertTrue[3=count pnlDist;"3 scenarios"];
pnlVec:.riskdist.pnlVector pnlDist;
.testutil.assertTrue[3=count pnlVec;"3 PnL values"];

lossVec:.riskdist.lossVector pnlVec;
.testutil.assertTrue[3=count lossVec;"3 loss values"];

summaryResult:.riskdist.distributionSummary pnlVec;
.testutil.assertTrue[summaryResult[`observationCount]=3;"3 observations"];
.testutil.assertTrue[not null summaryResult`meanPnl;"mean finite"];

worstRows:.riskdist.worstLossScenarios[pnlDist;2];
.testutil.assertTrue[2=count worstRows;"2 worst scenarios"];
worstRow:worstRows 0;
.testutil.assertTrue[worstRow[`rank]=1;"rank starts at 1"];
.testutil.assertTrue[worstRow[`loss]>=0f;"loss is positive"];

-1 "PASS test_pnl_distribution: scenarios=3, worstLoss=",string worstRow`loss;
