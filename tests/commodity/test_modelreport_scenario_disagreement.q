\l lib/init.q
/ scenarioDisagreement summarises the scenarioPnL table from scenarioShift and
/ raises scenarioPnlAlert when the PnL range across models exceeds threshold.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;10000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

modelInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

scenarioCfg:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);
basePriceTbl:.commodity.modelreport.priceAllModels[optionSetup;modelInputs];
scenPriceTbl:.commodity.modelreport.scenarioShift[optionSetup;modelInputs;scenarioCfg];
pnLTbl:.commodity.modelreport.scenarioPnL[basePriceTbl;scenPriceTbl;1f;1000f];

disCfg:.commodity.modelreport.defaultDisagreementConfig[];
scenDis:.commodity.modelreport.scenarioDisagreement[pnLTbl;disCfg];

expectedKeys:`okModelCount`minScenarioPnl`maxScenarioPnl`averageScenarioPnl`scenarioPnlRange`scenarioPnlAlert`scenarioPnlRangeThreshold`status`errorMessage;
.testutil.assertTrue[all expectedKeys in key scenDis;"scenarioDisagreement dict has expected keys"];
.testutil.assertTrue[4=scenDis`okModelCount;"okModelCount = 4"];
.testutil.assertTrue[`OK=scenDis`status;"status OK"];
.testutil.assertNear[scenDis`scenarioPnlRange;scenDis[`maxScenarioPnl]-scenDis`minScenarioPnl;1e-9;"PnL range = max - min"];
.testutil.assertTrue[-1h=type scenDis`scenarioPnlAlert;"alert boolean"];

expectedAlert:scenDis[`scenarioPnlRange]>disCfg`scenarioPnlRangeThreshold;
.testutil.assertTrue[scenDis[`scenarioPnlAlert]=expectedAlert;"alert flag matches threshold rule"];

strictCfg:@[disCfg;`scenarioPnlRangeThreshold;:;0.001];
strictDis:.commodity.modelreport.scenarioDisagreement[pnLTbl;strictCfg];
.testutil.assertTrue[strictDis`scenarioPnlAlert;"strict PnL threshold forces alert"];

-1 "PASS test_modelreport_scenario_disagreement: range=",string[scenDis`scenarioPnlRange],", alert=",string[scenDis`scenarioPnlAlert];
