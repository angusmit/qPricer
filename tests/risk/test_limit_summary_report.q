\l core/init.q
/ Build check result with mixed severities
checkResult:();
checkResult:checkResult,enlist `limitId`scopeType`scopeValue`metricName`metricValue`limitValue`warningThreshold`hardThreshold`utilisation`severity`passed`status`errorMessage!(1;`portfolio;`ALL;`VaR95;500000f;1000000f;800000f;1000000f;0.5;`OK;1b;`OK;"");
checkResult:checkResult,enlist `limitId`scopeType`scopeValue`metricName`metricValue`limitValue`warningThreshold`hardThreshold`utilisation`severity`passed`status`errorMessage!(2;`portfolio;`ALL;`ES95;900000f;1000000f;800000f;1000000f;0.9;`warning;1b;`OK;"");
checkResult:checkResult,enlist `limitId`scopeType`scopeValue`metricName`metricValue`limitValue`warningThreshold`hardThreshold`utilisation`severity`passed`status`errorMessage!(3;`book;`EQD;`DeltaCash;1200000f;1000000f;800000f;1000000f;1.2;`breach;0b;`OK;"");

breaches:.limitreport.breachReport checkResult;
.testutil.assertTrue[1=count breaches;"1 breach"];

warnings:.limitreport.warningReport checkResult;
.testutil.assertTrue[1=count warnings;"1 warning"];

sevSummary:.limitreport.summaryBySeverity checkResult;
.testutil.assertTrue[(count sevSummary)>0;"severity summary rows"];

scopeSummary:.limitreport.summaryByScope checkResult;
.testutil.assertTrue[(count scopeSummary)>0;"scope summary rows"];

metricSummary:.limitreport.summaryByMetric checkResult;
.testutil.assertTrue[(count metricSummary)>0;"metric summary rows"];

dashboard:.limitreport.limitDashboard checkResult;
.testutil.assertTrue[dashboard[`totalLimits]=3;"3 total limits"];
.testutil.assertTrue[dashboard[`okCount]=1;"1 OK"];
.testutil.assertTrue[dashboard[`warningCount]=1;"1 warning"];
.testutil.assertTrue[dashboard[`breachCount]=1;"1 breach"];
.testutil.assertTrue[dashboard[`worstSeverity]=`breach;"worst=breach"];

-1 "PASS test_limit_summary_report";
