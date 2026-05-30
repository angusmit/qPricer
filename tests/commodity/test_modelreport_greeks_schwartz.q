\l core/init.q
/ greeksSchwartz produces FD sensitivities for the Schwartz one-factor model.
/ Check schema, positivity of basePrice and volatilityVega, and that the
/ structural sensitivities are finite (mean-reversion and longRunLogMean
/ sensitivities are signed but their sign depends on the spot/theta relation,
/ so we only assert finiteness).

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);
modelInputs:enlist[`schwartz]!enlist `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));

greekCfg:.commodity.modelreport.defaultGreekConfig[];
greekTbl:.commodity.modelreport.greeksSchwartz[optionSetup;modelInputs;greekCfg];

requiredCols:`modelName`basePrice`logStateDelta`volatilityVega`meanReversionSensitivity`longRunLogMeanSensitivity`pricingStatus`pricingErrorMessage;
.testutil.assertTableColumns[greekTbl;requiredCols;"greeksSchwartz schema"];
.testutil.assertTrue[1=count greekTbl;"greeksSchwartz returns one row"];

greekRow:greekTbl 0;
.testutil.assertTrue[greekRow[`pricingStatus]=`OK;"schwartz greeks status OK"];
.testutil.assertTrue[greekRow[`basePrice]>0f;"basePrice positive"];
.testutil.assertTrue[not null greekRow`logStateDelta;"logStateDelta finite"];
.testutil.assertTrue[greekRow[`logStateDelta]>0f;"logStateDelta positive for call"];
.testutil.assertTrue[greekRow[`volatilityVega]>0f;"volatilityVega positive for call"];
.testutil.assertTrue[not null greekRow`meanReversionSensitivity;"meanReversionSensitivity finite"];
.testutil.assertTrue[not null greekRow`longRunLogMeanSensitivity;"longRunLogMeanSensitivity finite"];

-1 "PASS test_modelreport_greeks_schwartz: basePrice=",string[greekRow`basePrice],", delta=",string[greekRow`logStateDelta],", vega=",string[greekRow`volatilityVega];
