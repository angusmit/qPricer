\l lib/init.q
/ Synthetic reports
varReport:();
varReport:varReport,enlist `confidenceLevel`valueAtRisk`expectedShortfall`tailCount`observationCount`status`errorMessage!(0.95;400000f;500000f;5;100;`OK;"");
varReport:varReport,enlist `confidenceLevel`valueAtRisk`expectedShortfall`tailCount`observationCount`status`errorMessage!(0.99;700000f;800000f;1;100;`OK;"");

greekReport:enlist `scopeType`scopeValue`DeltaCash`VegaCash!(`portfolio;`ALL;250000f;180000f);
pnlReport:enlist `scenarioName`pnl`loss`rank`status`errorMessage!(`worst;-600000f;600000f;1;`OK;"");

metricTable:.limits.buildMetricTable[varReport;greekReport;pnlReport;(::)];
.testutil.assertTrue[(count metricTable)>0;"metric table has rows"];

limitTable:();
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`VaR95;500000f;0.8;1.0;`lessThan;1b);
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(2;`portfolio;`ALL;`VaR99;500000f;0.8;1.0;`lessThan;1b);
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(3;`portfolio;`ALL;`DeltaCash;300000f;0.8;1.0;`absLessThan;1b);
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(4;`portfolio;`ALL;`WorstLoss;500000f;0.8;1.0;`lessThan;1b);

checkResult:.limits.evaluateLimits[limitTable;metricTable];
.testutil.assertTrue[4=count checkResult;"4 limits checked"];

dashboard:.limitreport.limitDashboard checkResult;
.testutil.assertTrue[dashboard[`totalLimits]=4;"4 total"];
.testutil.assertTrue[dashboard[`okCount]>0;"some OK"];

-1 "PASS test_portfolio_limit_monitoring: total=",string[dashboard`totalLimits],", breaches=",string[dashboard`breachCount],", warnings=",string dashboard`warningCount;
