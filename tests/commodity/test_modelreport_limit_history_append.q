\l lib/init.q
/ appendLimitHistory accumulates snapshot rows. Verify schema stability,
/ row growth, runId preservation, and that an empty seed history returns
/ the new snapshot unchanged.

baseLimitTable:([] limitName:`grossPriceRangeExposure`maxScenarioPnlRange;
                   metricValue:4000 500f;
                   warningThreshold:5000 1000f;
                   breachThreshold:10000 2500f;
                   status:`OK`OK;
                   breachAmount:0 0f;
                   message:("";""));

snap1:.commodity.modelreport.limitSnapshot[`R1;2026.05.27;`bookA;baseLimitTable];
snap2:.commodity.modelreport.limitSnapshot[`R2;2026.05.28;`bookA;baseLimitTable];

seeded:.commodity.modelreport.appendLimitHistory[();snap1];
.testutil.assertTrue[(count seeded)=count snap1;"empty seed -> returns first snapshot"];
.testutil.assertTrue[(cols seeded)~cols snap1;"seeded schema matches snapshot schema"];

twoRunHistory:.commodity.modelreport.appendLimitHistory[seeded;snap2];
.testutil.assertTrue[(2*count snap1)=count twoRunHistory;"row count doubles after appending second snapshot"];
.testutil.assertTrue[(cols twoRunHistory)~cols snap1;"schema stable after append"];

uniqueRunIds:asc distinct twoRunHistory`runId;
.testutil.assertTrue[uniqueRunIds~asc `R1`R2;"both runIds present"];

snap3:.commodity.modelreport.limitSnapshot[`R3;2026.05.29;`bookB;baseLimitTable];
threeRunHistory:.commodity.modelreport.appendLimitHistory[twoRunHistory;snap3];
.testutil.assertTrue[(3*count snap1)=count threeRunHistory;"row count grows to 3 runs"];
.testutil.assertTrue[3=count distinct threeRunHistory`runId;"three distinct runIds"];
.testutil.assertTrue[2=count distinct threeRunHistory`portfolioName;"two distinct portfolios after appending bookB"];

emptyAppendedToHistory:.commodity.modelreport.appendLimitHistory[twoRunHistory;()];
.testutil.assertTrue[(count twoRunHistory)=count emptyAppendedToHistory;"empty new snapshot leaves history unchanged"];

-1 "PASS test_modelreport_limit_history_append: rowsAfter=",string[count threeRunHistory];
