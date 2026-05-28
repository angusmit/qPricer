\l lib/init.q
/ limitHistorySummary aggregates per-limitName statistics across runs.
/ Build a three-run history where one limit transitions OK -> warning ->
/ breach and the other stays OK. Verify counts, rates, latestStatus, and
/ maxBreachAmount.

mkSnap:{[runId;runDate;gMetric;sMetric]
    gStatus:?[gMetric>=10000;`breach;?[gMetric>=5000;`warning;`OK]];
    sStatus:?[sMetric>=2500;`breach;?[sMetric>=1000;`warning;`OK]];
    gBreach:?[gMetric>=10000;gMetric-10000;0f];
    sBreach:?[sMetric>=2500;sMetric-2500;0f];
    tbl:([] limitName:`grossPriceRangeExposure`maxScenarioPnlRange;
            metricValue:(`float$gMetric;`float$sMetric);
            warningThreshold:5000 1000f;
            breachThreshold:10000 2500f;
            status:(gStatus;sStatus);
            breachAmount:(gBreach;sBreach);
            message:("";""));
    .commodity.modelreport.limitSnapshot[runId;runDate;`bookA;tbl]
    };

run1:mkSnap[`R1;2026.05.27;4000f;500f];
run2:mkSnap[`R2;2026.05.28;6000f;500f];
run3:mkSnap[`R3;2026.05.29;12000f;500f];

history:.commodity.modelreport.appendLimitHistory[();run1];
history:.commodity.modelreport.appendLimitHistory[history;run2];
history:.commodity.modelreport.appendLimitHistory[history;run3];

historySum:.commodity.modelreport.limitHistorySummary history;

requiredCols:`limitName`observationCount`okCount`warningCount`breachCount`errorCount`latestStatus`latestMetricValue`worstMetricValue`maxBreachAmount`breachRate`warningRate`status`errorMessage;
.testutil.assertTableColumns[historySum;requiredCols;"limitHistorySummary schema"];
.testutil.assertTrue[2=count historySum;"two distinct limits in history"];

grossRow:historySum first where (historySum`limitName)=`grossPriceRangeExposure;
.testutil.assertTrue[3=grossRow`observationCount;"grossPriceRangeExposure observed 3 times"];
.testutil.assertTrue[1=grossRow`okCount;"grossPriceRangeExposure OK count = 1"];
.testutil.assertTrue[1=grossRow`warningCount;"grossPriceRangeExposure warning count = 1"];
.testutil.assertTrue[1=grossRow`breachCount;"grossPriceRangeExposure breach count = 1"];
.testutil.assertTrue[`breach=grossRow`latestStatus;"grossPriceRangeExposure latestStatus = breach (run3)"];
.testutil.assertNear[grossRow`latestMetricValue;12000f;1e-12;"grossPriceRangeExposure latest metric = 12000"];
.testutil.assertNear[grossRow`worstMetricValue;12000f;1e-12;"grossPriceRangeExposure worst = 12000"];
.testutil.assertNear[grossRow`maxBreachAmount;2000f;1e-12;"grossPriceRangeExposure maxBreachAmount = 12000 - 10000"];
.testutil.assertNear[grossRow`breachRate;1f%3f;1e-12;"grossPriceRangeExposure breachRate = 1/3"];
.testutil.assertNear[grossRow`warningRate;1f%3f;1e-12;"grossPriceRangeExposure warningRate = 1/3"];

scenarioRow:historySum first where (historySum`limitName)=`maxScenarioPnlRange;
.testutil.assertTrue[3=scenarioRow`okCount;"maxScenarioPnlRange OK count = 3"];
.testutil.assertTrue[0=scenarioRow`breachCount;"maxScenarioPnlRange breach count = 0"];
.testutil.assertTrue[`OK=scenarioRow`latestStatus;"maxScenarioPnlRange latestStatus = OK"];

emptySum:.commodity.modelreport.limitHistorySummary ();
.testutil.assertTrue[0=count emptySum;"empty history -> empty summary"];
.testutil.assertTableColumns[emptySum;requiredCols;"empty summary keeps schema"];

-1 "PASS test_modelreport_limit_history_summary: limits=",string[count historySum],", grossLatest=",string[grossRow`latestStatus];
