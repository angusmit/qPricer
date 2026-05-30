\l core/init.q
/ runComparisonRisk wraps runComparisonWithGreeks and modelDisagreementReport
/ into a single call. Verify the dictionary structure and that the alerts
/ table is populated with the expected alert types.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;8000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

modelInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

scenarioCfg:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);
greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

fullReport:.commodity.modelreport.runComparisonRisk[optionSetup;modelInputs;scenarioCfg;1f;1000f;greekCfg;disCfg];

.testutil.assertTrue[all `comparison`disagreement in key fullReport;"runComparisonRisk has comparison + disagreement"];

comparisonResult:fullReport`comparison;
.testutil.assertTrue[all `basePrices`greeksTable`scenarioPnL in key comparisonResult;"comparison has basePrices/greeksTable/scenarioPnL"];

disagreementResult:fullReport`disagreement;
disKeys:`priceDisagreement`greeksDisagreement`scenarioDisagreement`alerts;
.testutil.assertTrue[all disKeys in key disagreementResult;"disagreement has all four sub-keys"];

alertsTbl:disagreementResult`alerts;
.testutil.assertTrue[6=count alertsTbl;"alerts table has six rows"];
expectedAlertTypes:`priceRangeAbs`priceRangePct`primarySensitivityRange`volatilityVegaRange`scenarioPnlRange`jumpSensitivity;
.testutil.assertTrue[all expectedAlertTypes in alertsTbl`alertType;"all expected alert types present"];

priceDis:disagreementResult`priceDisagreement;
.testutil.assertTrue[`OK=priceDis`status;"priceDisagreement status OK"];
.testutil.assertTrue[priceDis[`okModelCount]>=2;"priceDisagreement has at least minimum OK models"];

greeksDis:disagreementResult`greeksDisagreement;
.testutil.assertTrue[`OK=greeksDis`status;"greeksDisagreement status OK"];

scenDis:disagreementResult`scenarioDisagreement;
.testutil.assertTrue[`OK=scenDis`status;"scenarioDisagreement status OK"];

-1 "PASS test_modelreport_comparison_risk: alertRows=",string[count alertsTbl],", warningCount=",string[sum (alertsTbl`severity)=`warning];
