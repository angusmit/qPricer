\l lib/init.q
/ Synthetic inputs
pricingResult:();
pricingResult:pricingResult,enlist `tradeId`underlying`optionType`unitPrice`notionalPrice`method`modelName`status!(1;`AAPL;`call;10.5;1050000f;`crankNicolson;`blackScholes;`OK);
pricingResult:pricingResult,enlist `tradeId`underlying`optionType`unitPrice`notionalPrice`method`modelName`status!(2;`MSFT;`call;20.0;1000000f;`crankNicolson;`blackScholes;`OK);

greekResult:();
greekResult:greekResult,enlist `delta`gamma`vega`theta`rho!(0.6;0.02;37.5;-6.4;53.2);
greekResult:greekResult,enlist `delta`gamma`vega`theta`rho!(0.5;0.01;25.0;-5.0;40.0);

varReport:enlist `confidenceLevel`valueAtRisk`expectedShortfall`tailCount`observationCount`status`errorMessage!(0.95;300000f;400000f;5;100;`OK;"");

limitDash:`totalLimits`okCount`warningCount`breachCount`errorCount`worstSeverity`maxUtilisation`status`errorMessage!(4;2;1;1;0;`breach;1.2;`OK;"");
mcSummary:`checkCount`passedCount`failedCount`maxAbsoluteDifference`maxRelativeDifference`status!(9;9;0;0.03;0.005;`OK);

dailyRiskResult:`pricingResult`greekResult`scenarioResult`pnlExplainResult`varReport`historicalReplayResult`historicalPnlDistribution`historicalVarReport`worstHistoricalEvents`limitCheckResult`limitDashboard`modelCheckResult`modelCheckSummary`runStatus!(
    pricingResult;greekResult;();();varReport;();();();();();limitDash;();mcSummary;());

dashSummary:.dashboard.dashboardSummary dailyRiskResult;
.testutil.assertTrue[`tradeCount in key dashSummary;"has tradeCount"];
.testutil.assertTrue[dashSummary[`tradeCount]=2;"2 trades"];
.testutil.assertTrue[not null dashSummary`totalPV;"totalPV populated"];
.testutil.assertTrue[not null dashSummary`var95;"var95 populated"];
.testutil.assertTrue[dashSummary[`limitBreachCount]=1;"1 breach"];
.testutil.assertTrue[dashSummary[`overallStatus]=`breach;"overall=breach"];
.testutil.assertTrue[dashSummary[`modelCheckPassedCount]=9;"9 model checks passed"];

-1 "PASS test_dashboard: totalPV=",string[dashSummary`totalPV],", overallStatus=",string dashSummary`overallStatus;
