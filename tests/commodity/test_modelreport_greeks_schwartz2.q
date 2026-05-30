\l core/init.q
/ greeksSchwartz2 produces FD sensitivities for the Schwartz two-factor model.
/ Long-factor moves persist through expiry (random walk), short-factor moves
/ are mean-reverting, so longFactorSensitivity is expected to be larger in
/ magnitude than shortFactorSensitivity for moderate kappa.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);
modelInputs:enlist[`schwartz2]!enlist `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));

greekCfg:.commodity.modelreport.defaultGreekConfig[];
greekTbl:.commodity.modelreport.greeksSchwartz2[optionSetup;modelInputs;greekCfg];

requiredCols:`modelName`basePrice`shortFactorSensitivity`longFactorSensitivity`shortVolSensitivity`longVolSensitivity`meanReversionSensitivity`correlationSensitivity`pricingStatus`pricingErrorMessage;
.testutil.assertTableColumns[greekTbl;requiredCols;"greeksSchwartz2 schema"];
.testutil.assertTrue[1=count greekTbl;"greeksSchwartz2 returns one row"];

greekRow:greekTbl 0;
.testutil.assertTrue[greekRow[`pricingStatus]=`OK;"schwartz2 greeks status OK"];
.testutil.assertTrue[greekRow[`basePrice]>0f;"basePrice positive"];
.testutil.assertTrue[not null greekRow`shortFactorSensitivity;"shortFactorSensitivity finite"];
.testutil.assertTrue[not null greekRow`longFactorSensitivity;"longFactorSensitivity finite"];
.testutil.assertTrue[not null greekRow`shortVolSensitivity;"shortVolSensitivity finite"];
.testutil.assertTrue[not null greekRow`longVolSensitivity;"longVolSensitivity finite"];
.testutil.assertTrue[not null greekRow`meanReversionSensitivity;"meanReversionSensitivity finite"];
.testutil.assertTrue[not null greekRow`correlationSensitivity;"correlationSensitivity finite"];
.testutil.assertTrue[greekRow[`shortFactorSensitivity]>0f;"shortFactorSensitivity positive for call"];
.testutil.assertTrue[greekRow[`longFactorSensitivity]>0f;"longFactorSensitivity positive for call"];
.testutil.assertTrue[greekRow[`longVolSensitivity]>0f;"longVolSensitivity positive for call"];
.testutil.assertTrue[greekRow[`longFactorSensitivity]>greekRow`shortFactorSensitivity;"long factor sensitivity dominates short factor (persists past mean reversion)"];

-1 "PASS test_modelreport_greeks_schwartz2: short=",string[greekRow`shortFactorSensitivity],", long=",string[greekRow`longFactorSensitivity],", shortVol=",string[greekRow`shortVolSensitivity],", longVol=",string[greekRow`longVolSensitivity];
