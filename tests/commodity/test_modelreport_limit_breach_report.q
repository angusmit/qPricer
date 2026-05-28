\l lib/init.q
/ limitBreachReport excludes OK rows and sorts by severity (ERROR, breach,
/ warning) then by breachAmount descending. Verify ordering on a synthetic
/ limit table that mixes all four statuses, including a non-breach OK row
/ that must be excluded.

okRow:`limitName`metricValue`warningThreshold`breachThreshold`status`breachAmount`message!(
    `okOne;100f;500f;1000f;`OK;0f;"");
warnRow:`limitName`metricValue`warningThreshold`breachThreshold`status`breachAmount`message!(
    `warnOne;600f;500f;1000f;`warning;0f;"warn");
breachSmall:`limitName`metricValue`warningThreshold`breachThreshold`status`breachAmount`message!(
    `breachSmall;1100f;500f;1000f;`breach;100f;"small");
breachLarge:`limitName`metricValue`warningThreshold`breachThreshold`status`breachAmount`message!(
    `breachLarge;1500f;500f;1000f;`breach;500f;"large");
errorRow:`limitName`metricValue`warningThreshold`breachThreshold`status`breachAmount`message!(
    `errOne;0Nf;500f;1000f;`ERROR;0Nf;"err");

synthTable:(,/) enlist each (okRow;warnRow;breachSmall;breachLarge;errorRow);
report:.commodity.modelreport.limitBreachReport synthTable;

requiredCols:`limitName`metricValue`warningThreshold`breachThreshold`status`breachAmount`message`breachRank;
.testutil.assertTableColumns[report;requiredCols;"limitBreachReport schema includes breachRank"];

.testutil.assertTrue[4=count report;"OK row excluded (5 -> 4)"];
.testutil.assertTrue[not `okOne in report`limitName;"okOne explicitly absent"];

reportStatuses:report`status;
.testutil.assertTrue[`ERROR=reportStatuses 0;"first row is ERROR"];
.testutil.assertTrue[`breach=reportStatuses 1;"second row is breach"];
.testutil.assertTrue[`breach=reportStatuses 2;"third row is breach"];
.testutil.assertTrue[`warning=reportStatuses 3;"fourth row is warning"];

.testutil.assertTrue[`breachLarge=report[1;`limitName];"larger breachAmount ranked before smaller within breach group"];
.testutil.assertTrue[`breachSmall=report[2;`limitName];"smaller breachAmount second within breach group"];

ranksCol:report`breachRank;
.testutil.assertTrue[ranksCol~1+til 4;"breachRank = 1..N"];

emptyReport:.commodity.modelreport.limitBreachReport ();
.testutil.assertTableColumns[emptyReport;requiredCols;"empty input keeps schema"];
.testutil.assertTrue[0=count emptyReport;"empty input -> empty report"];

allOkTable:(,/) enlist each (okRow;@[okRow;`limitName;:;`okTwo]);
allOkReport:.commodity.modelreport.limitBreachReport allOkTable;
.testutil.assertTrue[0=count allOkReport;"all-OK input -> empty report"];

-1 "PASS test_modelreport_limit_breach_report: rows=",string[count report],", topStatus=",string reportStatuses 0;
