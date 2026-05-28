\l lib/init.q
/ disagreementAlerts builds a six-row table from already-computed disagreement
/ dicts. Feed it controlled synthetic dicts so the test does not depend on the
/ pricing models. Verify alert types, severity mapping, and the metric/threshold
/ columns.

priceDisAlert:`okModelCount`minPrice`maxPrice`averagePrice`priceRange`priceRangePct`priceRangeAlert`priceRangeAbsThreshold`priceRangePctThreshold`status`errorMessage!(
    4;1f;10f;5f;9f;1.8;1b;5f;0.10;`OK;"");
priceDisQuiet:`okModelCount`minPrice`maxPrice`averagePrice`priceRange`priceRangePct`priceRangeAlert`priceRangeAbsThreshold`priceRangePctThreshold`status`errorMessage!(
    4;9.99f;10.01f;10f;0.02f;0.002;0b;5f;0.10;`OK;"");

greeksDis:`okModelCount`primarySensitivityMin`primarySensitivityMax`primarySensitivityRange`primarySensitivityAlert`volatilityVegaMin`volatilityVegaMax`volatilityVegaRange`volatilityVegaAlert`jumpIntensitySensitivity`jumpSensitivityAlert`primaryDeltaRangeThreshold`volatilityVegaRangeThreshold`jumpSensitivityThreshold`status`errorMessage!(
    4;0.5;44f;43.5;1b;15f;28f;13f;1b;1.9;0b;10f;10f;5f;`OK;"");

scenDis:`okModelCount`minScenarioPnl`maxScenarioPnl`averageScenarioPnl`scenarioPnlRange`scenarioPnlAlert`scenarioPnlRangeThreshold`status`errorMessage!(
    4;1000f;3500f;2200f;2500f;1b;1000f;`OK;"");

alertsAlert:.commodity.modelreport.disagreementAlerts[priceDisAlert;greeksDis;scenDis];

requiredCols:`alertType`alertFlag`metricValue`threshold`severity`message;
.testutil.assertTableColumns[alertsAlert;requiredCols;"disagreementAlerts schema"];

expectedAlertTypes:`priceRangeAbs`priceRangePct`primarySensitivityRange`volatilityVegaRange`scenarioPnlRange`jumpSensitivity;
.testutil.assertTrue[all expectedAlertTypes in alertsAlert`alertType;"all six alert types present"];
.testutil.assertTrue[6=count alertsAlert;"six alert rows"];

severityByType:exec severity by alertType from alertsAlert;
.testutil.assertTrue[`warning=first severityByType`priceRangeAbs;"priceRangeAbs severity warning (range 9 > threshold 5)"];
.testutil.assertTrue[`warning=first severityByType`priceRangePct;"priceRangePct severity warning"];
.testutil.assertTrue[`warning=first severityByType`primarySensitivityRange;"primarySensitivityRange severity warning"];
.testutil.assertTrue[`warning=first severityByType`volatilityVegaRange;"volatilityVegaRange severity warning"];
.testutil.assertTrue[`warning=first severityByType`scenarioPnlRange;"scenarioPnlRange severity warning"];
.testutil.assertTrue[`OK=first severityByType`jumpSensitivity;"jumpSensitivity severity OK (1.9 < threshold 5)"];

alertsQuiet:.commodity.modelreport.disagreementAlerts[priceDisQuiet;greeksDis;scenDis];
quietByType:exec severity by alertType from alertsQuiet;
.testutil.assertTrue[`OK=first quietByType`priceRangeAbs;"quiet price abs severity OK"];
.testutil.assertTrue[`OK=first quietByType`priceRangePct;"quiet price pct severity OK"];

allBool:all -1h=type each alertsAlert`alertFlag;
.testutil.assertTrue[allBool;"alertFlag column is boolean"];

-1 "PASS test_modelreport_disagreement_alerts: rows=",string[count alertsAlert],", warningCount=",string[sum (alertsAlert`severity)=`warning];
