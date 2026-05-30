\l core/init.q
/ greeksDisagreement summarises greeksAllModels into primary-sensitivity range,
/ volatility-vega range, and jump-intensity sensitivity for mrjump. Each metric
/ has a configurable threshold; alerts must be boolean.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;8000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

modelInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

greekCfg:.commodity.modelreport.defaultGreekConfig[];
greeksTbl:.commodity.modelreport.greeksAllModels[optionSetup;modelInputs;greekCfg];

disCfg:.commodity.modelreport.defaultDisagreementConfig[];
greeksDis:.commodity.modelreport.greeksDisagreement[greeksTbl;disCfg];

expectedKeys:`okModelCount`primarySensitivityMin`primarySensitivityMax`primarySensitivityRange`primarySensitivityAlert`volatilityVegaMin`volatilityVegaMax`volatilityVegaRange`volatilityVegaAlert`jumpIntensitySensitivity`jumpSensitivityAlert`primaryDeltaRangeThreshold`volatilityVegaRangeThreshold`jumpSensitivityThreshold`status`errorMessage;
.testutil.assertTrue[all expectedKeys in key greeksDis;"greeksDisagreement dict has expected keys"];

.testutil.assertTrue[greeksDis[`primarySensitivityRange]>=0f;"primarySensitivityRange non-negative"];
.testutil.assertTrue[greeksDis[`volatilityVegaRange]>=0f;"volatilityVegaRange non-negative"];
.testutil.assertNear[greeksDis`primarySensitivityRange;greeksDis[`primarySensitivityMax]-greeksDis`primarySensitivityMin;1e-12;"primary range = max - min"];
.testutil.assertNear[greeksDis`volatilityVegaRange;greeksDis[`volatilityVegaMax]-greeksDis`volatilityVegaMin;1e-12;"vega range = max - min"];

.testutil.assertTrue[-1h=type greeksDis`primarySensitivityAlert;"primary alert boolean"];
.testutil.assertTrue[-1h=type greeksDis`volatilityVegaAlert;"vega alert boolean"];
.testutil.assertTrue[-1h=type greeksDis`jumpSensitivityAlert;"jump alert boolean"];

.testutil.assertTrue[not null greeksDis`jumpIntensitySensitivity;"jumpIntensitySensitivity populated from mrjump row"];

greeksTblExMrJump:greeksTbl where not (greeksTbl`modelName)=`mrjump;
disExMrJump:.commodity.modelreport.greeksDisagreement[greeksTblExMrJump;disCfg];
.testutil.assertTrue[null disExMrJump`jumpIntensitySensitivity;"jumpIntensitySensitivity null when mrjump absent"];
.testutil.assertTrue[not disExMrJump`jumpSensitivityAlert;"jump alert false when mrjump absent"];

strictCfg:@[disCfg;`primaryDeltaRangeThreshold;:;0.001];
strictDis:.commodity.modelreport.greeksDisagreement[greeksTbl;strictCfg];
.testutil.assertTrue[strictDis`primarySensitivityAlert;"strict primary threshold forces alert"];

-1 "PASS test_modelreport_greeks_disagreement: primaryRange=",string[greeksDis`primarySensitivityRange],", vegaRange=",string[greeksDis`volatilityVegaRange],", jumpSens=",string[greeksDis`jumpIntensitySensitivity];
