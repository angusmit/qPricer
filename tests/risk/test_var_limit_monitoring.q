\l core/init.q
/ Synthetic VaR report
varReport:();
varReport:varReport,enlist `confidenceLevel`valueAtRisk`expectedShortfall`tailCount`observationCount`status`errorMessage!(0.95;500000f;600000f;5;100;`OK;"");

/ Limits: VaR95 limit = 1M (OK), ES95 limit = 400K (breach)
limitTable:();
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`VaR95;1000000f;0.8;1.0;`lessThan;1b);
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(2;`portfolio;`ALL;`ES95;400000f;0.8;1.0;`lessThan;1b);

checkResult:.limits.checkVarLimits[varReport;limitTable];
.testutil.assertTrue[2=count checkResult;"2 limit checks"];

/ VaR95=500K vs limit=1M: utilisation=0.5 -> OK
varRow:checkResult 0;
.testutil.assertTrue[varRow[`severity]=`OK;"VaR95 OK"];
.testutil.assertNear[varRow`utilisation;0.5;0.01;"VaR utilisation ~0.5"];

/ ES95=600K vs limit=400K: 600K > 400K -> breach
esRow:checkResult 1;
.testutil.assertTrue[esRow[`severity]=`breach;"ES95 breach"];

-1 "PASS test_var_limit_monitoring";
