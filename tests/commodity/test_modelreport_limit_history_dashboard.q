\l lib/init.q
/ limitHistoryDashboard composes historySummary + breachTrend +
/ repeatedBreaches + latestSummary into a single dashboard dict.

mkSnap:{[runId;runDate;gMetric;status]
    breach:?[status=`breach;gMetric-10000;0f];
    tbl:([] limitName:enlist `grossPriceRangeExposure;
            metricValue:enlist `float$gMetric;
            warningThreshold:enlist 5000f;
            breachThreshold:enlist 10000f;
            status:enlist status;
            breachAmount:enlist breach;
            message:enlist "");
    .commodity.modelreport.limitSnapshot[runId;runDate;`bookD;tbl]
    };

mkSumSnap:{[runId;runDate;overall]
    summaryDict:`limitCount`okCount`warningCount`breachCount`errorCount`overallStatus`statusMessage!(
        1;?[overall=`OK;1;0];?[overall=`warning;1;0];?[overall=`breach;1;0];0;overall;"x");
    .commodity.modelreport.limitSummarySnapshot[runId;runDate;`bookD;summaryDict]
    };

run1:mkSnap[`R1;2026.05.27;4000f;`OK];
run2:mkSnap[`R2;2026.05.28;7000f;`warning];
run3:mkSnap[`R3;2026.05.29;12000f;`breach];
sum1:mkSumSnap[`R1;2026.05.27;`OK];
sum2:mkSumSnap[`R2;2026.05.28;`warning];
sum3:mkSumSnap[`R3;2026.05.29;`breach];

limitHistory:.commodity.modelreport.appendLimitHistory[();run1];
limitHistory:.commodity.modelreport.appendLimitHistory[limitHistory;run2];
limitHistory:.commodity.modelreport.appendLimitHistory[limitHistory;run3];
summaryHistory:.commodity.modelreport.appendLimitHistory[();sum1];
summaryHistory:.commodity.modelreport.appendLimitHistory[summaryHistory;sum2];
summaryHistory:.commodity.modelreport.appendLimitHistory[summaryHistory;sum3];

dashboard:.commodity.modelreport.limitHistoryDashboard[limitHistory;summaryHistory;2;1];

expectedKeys:`historySummary`breachTrend`repeatedBreaches`latestSummary`dashboardStatus`dashboardMessage;
.testutil.assertTrue[all expectedKeys in key dashboard;"dashboard has expected keys"];

historySum:dashboard`historySummary;
.testutil.assertTrue[1=count historySum;"one limit in historySummary"];

trendTbl:dashboard`breachTrend;
.testutil.assertTrue[1=count trendTbl;"one limit in trend"];
.testutil.assertTrue[`worsening=first trendTbl`trendDirection;"trend is worsening across the three runs"];

repeated:dashboard`repeatedBreaches;
.testutil.assertTrue[1=count repeated;"repeated breach surfaced with minBreachCount=1"];

latestSum:dashboard`latestSummary;
.testutil.assertTrue[99h=type latestSum;"latestSummary is a dict (one row of summaryHistory)"];
.testutil.assertTrue[`R3=latestSum`runId;"latestSummary picks the latest runId"];
.testutil.assertTrue[`breach=latestSum`overallStatus;"latestSummary picks the latest overallStatus"];

.testutil.assertTrue[`breach=dashboard`dashboardStatus;"dashboardStatus mirrors latest overallStatus"];

emptyDashboard:.commodity.modelreport.limitHistoryDashboard[();();2;1];
.testutil.assertTrue[`ERROR=emptyDashboard`dashboardStatus;"empty history -> dashboardStatus ERROR"];

-1 "PASS test_modelreport_limit_history_dashboard: dashStatus=",string[dashboard`dashboardStatus],", trendDir=",string[first trendTbl`trendDirection];
