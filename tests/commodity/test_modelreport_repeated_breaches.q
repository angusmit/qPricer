\l core/init.q
/ repeatedBreaches surfaces limits with breachCount >= minBreachCount.

mkSnap:{[runId;runDate;gMetric]
    gStatus:?[gMetric>=10000;`breach;?[gMetric>=5000;`warning;`OK]];
    gBreach:?[gMetric>=10000;gMetric-10000;0f];
    tbl:([] limitName:enlist `grossPriceRangeExposure;
            metricValue:enlist `float$gMetric;
            warningThreshold:enlist 5000f;
            breachThreshold:enlist 10000f;
            status:enlist gStatus;
            breachAmount:enlist gBreach;
            message:enlist "");
    .commodity.modelreport.limitSnapshot[runId;runDate;`bookA;tbl]
    };

run1:mkSnap[`R1;2026.05.27;12000f];
run2:mkSnap[`R2;2026.05.28;11000f];
run3:mkSnap[`R3;2026.05.29;15000f];
run4:mkSnap[`R4;2026.05.30;4000f];

history:.commodity.modelreport.appendLimitHistory[();run1];
history:.commodity.modelreport.appendLimitHistory[history;run2];
history:.commodity.modelreport.appendLimitHistory[history;run3];
history:.commodity.modelreport.appendLimitHistory[history;run4];

repeated:.commodity.modelreport.repeatedBreaches[history;2];

requiredCols:`limitName`breachCount`warningCount`observationCount`breachRate`latestStatus`maxBreachAmount`status`errorMessage;
.testutil.assertTableColumns[repeated;requiredCols;"repeatedBreaches schema"];
.testutil.assertTrue[1=count repeated;"one limit qualifies as repeated breach"];
.testutil.assertTrue[`grossPriceRangeExposure=first repeated`limitName;"correct limit returned"];
.testutil.assertTrue[3=first repeated`breachCount;"breachCount = 3"];
.testutil.assertNear[first repeated`maxBreachAmount;5000f;1e-12;"maxBreachAmount = 15000 - 10000"];

repeatedHigh:.commodity.modelreport.repeatedBreaches[history;5];
.testutil.assertTrue[0=count repeatedHigh;"minBreachCount too high -> empty"];
.testutil.assertTableColumns[repeatedHigh;requiredCols;"empty repeated keeps schema"];

emptyRepeated:.commodity.modelreport.repeatedBreaches[();1];
.testutil.assertTrue[0=count emptyRepeated;"empty history -> empty repeated"];
.testutil.assertTableColumns[emptyRepeated;requiredCols;"empty history keeps schema"];

-1 "PASS test_modelreport_repeated_breaches: count=",string[count repeated],", breachCount=",string first repeated`breachCount;
