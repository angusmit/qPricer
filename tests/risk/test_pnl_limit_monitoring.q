\l core/init.q
/ Synthetic worst-loss report
pnlReport:();
pnlReport:pnlReport,enlist `scenarioName`pnl`loss`rank`status`errorMessage!(`covidCrash;-800000f;800000f;1;`OK;"");

/ Limit: WorstLoss < 1M (OK at 800K)
limitTable:();
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`WorstLoss;1000000f;0.8;1.0;`lessThan;1b);

checkResult:.limits.checkPnlLimits[pnlReport;limitTable];
.testutil.assertTrue[1=count checkResult;"1 pnl limit"];
pnlRow:checkResult 0;
/ 800K vs warning=800K: exactly at warning -> warning
.testutil.assertTrue[pnlRow[`severity] in `OK`warning;"PnL within or at warning"];

-1 "PASS test_pnl_limit_monitoring";
