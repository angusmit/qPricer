\l core/init.q
/ Synthetic Greek report with scope
greekReport:();
greekReport:greekReport,enlist `scopeType`scopeValue`DeltaCash`VegaCash!(`book;`EQD;-300000f;150000f);
greekReport:greekReport,enlist `scopeType`scopeValue`DeltaCash`VegaCash!(`underlying;`AAPL;120000f;80000f);

/ Limits: absLessThan
limitTable:();
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`book;`EQD;`DeltaCash;500000f;0.8;1.0;`absLessThan;1b);
limitTable:limitTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(2;`book;`EQD;`VegaCash;100000f;0.8;1.0;`absLessThan;1b);

checkResult:.limits.checkGreekLimits[greekReport;limitTable];
.testutil.assertTrue[2=count checkResult;"2 greek limits"];

/ DeltaCash=-300K, abs=300K vs limit 500K: OK
deltaRow:checkResult 0;
.testutil.assertTrue[deltaRow[`severity]=`OK;"delta OK"];

/ VegaCash=150K, abs=150K vs limit 100K: breach
vegaRow:checkResult 1;
.testutil.assertTrue[vegaRow[`severity]=`breach;"vega breach"];

-1 "PASS test_greek_limit_monitoring";
