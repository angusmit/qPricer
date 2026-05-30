\l core/init.q
scenarioResult:();
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`s1;-10f);
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`s2;-5f);
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`s3;2f);
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`s4;8f);
scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(`s5;-3f);

pnlDist:.riskdist.buildPnlDistribution scenarioResult;
worstRows:.riskdist.worstLossScenarios[pnlDist;3];
.testutil.assertTrue[3=count worstRows;"top 3 returned"];
.testutil.assertTrue[(worstRows 0)[`rank]=1;"rank 1"];
.testutil.assertTrue[(worstRows 0)[`loss]>=(worstRows 1)`loss;"sorted by loss"];
.testutil.assertTrue[(worstRows 0)[`loss]>0f;"loss positive"];

/ Requesting more than available
allRows:.riskdist.worstLossScenarios[pnlDist;100];
.testutil.assertTrue[5=count allRows;"capped at 5"];

-1 "PASS test_worst_loss_report: worstLoss=",string (worstRows 0)`loss;
