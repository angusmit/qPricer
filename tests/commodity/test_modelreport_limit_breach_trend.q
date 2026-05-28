\l lib/init.q
/ limitBreachTrend compares latest metric vs average of prior lookback runs
/ per limit. Verify worsening / improving / flat detection plus the
/ insufficient-history warning and the invalid-lookback error.

mkSnap:{[runId;runDate;gMetric;sMetric]
    tbl:([] limitName:`limitUp`limitDown;
            metricValue:(`float$gMetric;`float$sMetric);
            warningThreshold:5000 5000f;
            breachThreshold:10000 10000f;
            status:`OK`OK;
            breachAmount:0 0f;
            message:("";""));
    .commodity.modelreport.limitSnapshot[runId;runDate;`bookA;tbl]
    };

run1:mkSnap[`R1;2026.05.27;1000f;5000f];
run2:mkSnap[`R2;2026.05.28;2000f;4000f];
run3:mkSnap[`R3;2026.05.29;4000f;3000f];

history:.commodity.modelreport.appendLimitHistory[();run1];
history:.commodity.modelreport.appendLimitHistory[history;run2];
history:.commodity.modelreport.appendLimitHistory[history;run3];

trendTbl:.commodity.modelreport.limitBreachTrend[history;2];

requiredCols:`limitName`lookbackRuns`latestMetricValue`previousMetricValue`metricChange`metricChangePct`latestStatus`trendDirection`status`errorMessage;
.testutil.assertTableColumns[trendTbl;requiredCols;"limitBreachTrend schema"];
.testutil.assertTrue[2=count trendTbl;"two limits in trend table"];

upRow:trendTbl first where (trendTbl`limitName)=`limitUp;
.testutil.assertNear[upRow`latestMetricValue;4000f;1e-12;"limitUp latest = 4000"];
.testutil.assertNear[upRow`previousMetricValue;1500f;1e-12;"limitUp previous = avg(1000,2000)"];
.testutil.assertNear[upRow`metricChange;2500f;1e-12;"limitUp change = 2500"];
.testutil.assertTrue[`worsening=upRow`trendDirection;"limitUp trend = worsening"];

downRow:trendTbl first where (trendTbl`limitName)=`limitDown;
.testutil.assertNear[downRow`previousMetricValue;4500f;1e-12;"limitDown previous = avg(5000,4000)"];
.testutil.assertTrue[`improving=downRow`trendDirection;"limitDown trend = improving"];

flatSnap:mkSnap[`R4;2026.05.30;2000f;3000f];
flatHistory:.commodity.modelreport.appendLimitHistory[history;flatSnap];
flatHistory:.commodity.modelreport.appendLimitHistory[flatHistory;flatSnap];
flatHistoryTwoOnly:flatHistory where (flatHistory`runId) in `R4`R5;
flatTrend:.commodity.modelreport.limitBreachTrend[flatHistoryTwoOnly;1];
upFlatRow:flatTrend first where (flatTrend`limitName)=`limitUp;
.testutil.assertTrue[`flat=upFlatRow`trendDirection;"two identical runs -> flat"];

singleRun:.commodity.modelreport.appendLimitHistory[();run1];
singleTrend:.commodity.modelreport.limitBreachTrend[singleRun;2];
singleUpRow:singleTrend first where (singleTrend`limitName)=`limitUp;
.testutil.assertTrue[`flat=singleUpRow`trendDirection;"single-run trend = flat"];
.testutil.assertTrue[`warning=singleUpRow`status;"single-run trend status = warning"];

badLookback:@[.commodity.modelreport.limitBreachTrend[history;];0;{`ERROR}];
.testutil.assertTrue[badLookback~`ERROR;"lookback < 1 rejected"];

emptyTrend:.commodity.modelreport.limitBreachTrend[();2];
.testutil.assertTrue[0=count emptyTrend;"empty history -> empty trend"];

-1 "PASS test_modelreport_limit_breach_trend: upChange=",string[upRow`metricChange],", downDir=",string[downRow`trendDirection];
