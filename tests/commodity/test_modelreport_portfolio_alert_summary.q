\l core/init.q
/ portfolioAlertSummary aggregates alertRows by alertType. Test both a
/ synthetic input that controls counts exactly and an empty-input edge case.

okRow:`alertType`alertFlag`metricValue`threshold`severity`message`tradeId`commodity!(
    `priceRangeAbs;0b;3.0;5.0;`OK;"x";`T1;`wti);
warnRow:`alertType`alertFlag`metricValue`threshold`severity`message`tradeId`commodity!(
    `priceRangeAbs;1b;9.0;5.0;`warning;"x";`T2;`brent);
warnRow2:@[warnRow;(`alertType;`tradeId);:;(`scenarioPnlRange;`T3)];
synthAlertRows:(,/) enlist each (okRow;warnRow;warnRow2;@[okRow;`tradeId;:;`T4]);

alertSummary:.commodity.modelreport.portfolioAlertSummary synthAlertRows;

requiredCols:`alertType`okCount`warningCount`errorCount`totalCount`warningRate;
.testutil.assertTableColumns[alertSummary;requiredCols;"alertSummary schema"];
.testutil.assertTrue[2=count alertSummary;"two distinct alertTypes"];

priceAbsRow:alertSummary first where (alertSummary`alertType)=`priceRangeAbs;
.testutil.assertTrue[3=priceAbsRow`totalCount;"priceRangeAbs total = 3 rows"];
.testutil.assertTrue[2=priceAbsRow`okCount;"priceRangeAbs OK count = 2"];
.testutil.assertTrue[1=priceAbsRow`warningCount;"priceRangeAbs warning count = 1"];
.testutil.assertTrue[0=priceAbsRow`errorCount;"priceRangeAbs error count = 0"];
.testutil.assertNear[priceAbsRow`warningRate;1f%3f;1e-12;"priceRangeAbs warningRate = 1/3"];

scenarioRow:alertSummary first where (alertSummary`alertType)=`scenarioPnlRange;
.testutil.assertTrue[1=scenarioRow`totalCount;"scenarioPnlRange total = 1"];
.testutil.assertNear[scenarioRow`warningRate;1f;1e-12;"scenarioPnlRange warningRate = 1"];

emptySummary:.commodity.modelreport.portfolioAlertSummary ();
.testutil.assertTrue[0=count emptySummary;"empty alertRows -> empty summary"];
.testutil.assertTableColumns[emptySummary;requiredCols;"empty summary keeps schema"];

-1 "PASS test_modelreport_portfolio_alert_summary: rows=",string[count alertSummary];
