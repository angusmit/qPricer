\l core/init.q
/ limitSnapshot tags a limitTable with runId / runDate / portfolioName and
/ preserves row count. limitSummarySnapshot wraps a limitSummary dict into a
/ one-row table with the same metadata.

baseLimitTable:([] limitName:`grossPriceRangeExposure`maxScenarioPnlRange;
                   metricValue:4000 500f;
                   warningThreshold:5000 1000f;
                   breachThreshold:10000 2500f;
                   status:`OK`OK;
                   breachAmount:0 0f;
                   message:("";""));

snapshot:.commodity.modelreport.limitSnapshot[`R1;2026.05.27;`bookA;baseLimitTable];

requiredCols:`runId`runDate`portfolioName`limitName`metricValue`warningThreshold`breachThreshold`status`breachAmount`message;
.testutil.assertTableColumns[snapshot;requiredCols;"limitSnapshot schema"];
.testutil.assertTrue[2=count snapshot;"limitSnapshot row count preserved"];
.testutil.assertTrue[all (snapshot`runId)=`R1;"runId tag broadcast"];
.testutil.assertTrue[all (snapshot`runDate)=2026.05.27;"runDate tag broadcast"];
.testutil.assertTrue[all (snapshot`portfolioName)=`bookA;"portfolioName tag broadcast"];
.testutil.assertTrue[`runId=first cols snapshot;"runId column listed first"];

emptySnapshot:.commodity.modelreport.limitSnapshot[`R1;2026.05.27;`bookA;()];
.testutil.assertTrue[0=count emptySnapshot;"empty limitTable -> empty snapshot"];
.testutil.assertTableColumns[emptySnapshot;requiredCols;"empty snapshot keeps schema"];

baseSummary:`limitCount`okCount`warningCount`breachCount`errorCount`overallStatus`statusMessage!(7;5;1;1;0;`breach;"1 limit breach");
summarySnap:.commodity.modelreport.limitSummarySnapshot[`R1;2026.05.27;`bookA;baseSummary];

summaryCols:`runId`runDate`portfolioName`limitCount`okCount`warningCount`breachCount`errorCount`overallStatus`statusMessage;
.testutil.assertTableColumns[summarySnap;summaryCols;"summary snapshot schema"];
.testutil.assertTrue[1=count summarySnap;"summary snapshot one row"];
.testutil.assertTrue[`R1=first summarySnap`runId;"summary runId tag"];
.testutil.assertTrue[`breach=first summarySnap`overallStatus;"overallStatus passed through"];

-1 "PASS test_modelreport_limit_snapshot: snapshotRows=",string[count snapshot],", summaryRows=",string[count summarySnap];
