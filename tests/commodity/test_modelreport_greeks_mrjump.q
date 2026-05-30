\l core/init.q
/ greeksMrJump produces FD sensitivities for the mean-reverting jump model
/ using common-random-numbers (same mcConfig randomSeed) across base/up/down.
/ With moderate paths the structural checks (finite, positive standardError,
/ positive call vega) are stable; sign of mean-reversion and longRunLogMean
/ sensitivities depends on parameter regime so only finiteness is asserted.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;10000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

modelInputs:enlist[`mrjump]!enlist `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg);

greekCfg:.commodity.modelreport.defaultGreekConfig[];
greekTbl:.commodity.modelreport.greeksMrJump[optionSetup;modelInputs;greekCfg];

requiredCols:`modelName`basePrice`logStateDelta`volatilityVega`meanReversionSensitivity`longRunLogMeanSensitivity`jumpIntensitySensitivity`jumpMeanSensitivity`standardError`pricingStatus`pricingErrorMessage;
.testutil.assertTableColumns[greekTbl;requiredCols;"greeksMrJump schema"];
.testutil.assertTrue[1=count greekTbl;"greeksMrJump returns one row"];

greekRow:greekTbl 0;
.testutil.assertTrue[greekRow[`pricingStatus]=`OK;"mrjump greeks status OK"];
.testutil.assertTrue[greekRow[`basePrice]>0f;"basePrice positive"];
.testutil.assertTrue[greekRow[`standardError]>0f;"standardError positive"];
.testutil.assertTrue[not null greekRow`logStateDelta;"logStateDelta finite"];
.testutil.assertTrue[greekRow[`logStateDelta]>0f;"logStateDelta positive for call (spot up lifts price)"];
.testutil.assertTrue[not null greekRow`volatilityVega;"volatilityVega finite"];
.testutil.assertTrue[greekRow[`volatilityVega]>0f;"volatilityVega positive for call"];
.testutil.assertTrue[not null greekRow`meanReversionSensitivity;"meanReversionSensitivity finite"];
.testutil.assertTrue[not null greekRow`longRunLogMeanSensitivity;"longRunLogMeanSensitivity finite"];
.testutil.assertTrue[not null greekRow`jumpIntensitySensitivity;"jumpIntensitySensitivity finite"];
.testutil.assertTrue[not null greekRow`jumpMeanSensitivity;"jumpMeanSensitivity finite"];
.testutil.assertTrue[greekRow[`jumpMeanSensitivity]>0f;"jumpMeanSensitivity positive for call (positive jumps lift price)"];

-1 "PASS test_modelreport_greeks_mrjump: basePrice=",string[greekRow`basePrice]," SE=",string[greekRow`standardError]," delta=",string[greekRow`logStateDelta]," vega=",string[greekRow`volatilityVega]," jumpIntenSens=",string[greekRow`jumpIntensitySensitivity]," jumpMeanSens=",string[greekRow`jumpMeanSensitivity];
