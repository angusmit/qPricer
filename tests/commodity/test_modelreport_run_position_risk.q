\l core/init.q
/ runPositionRisk runs runComparisonRisk for one position and tags the
/ resulting tables with tradeId/commodity. Verify the returned dict shape and
/ that every tagged row carries the expected tradeId.

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;3000];
mcCfg:@[mcCfg;`timeStepCount;:;15];

setupCall:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);
inputsAll:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));
scenCfg:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);

position:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`T_PR;`wti;setupCall;inputsAll;scenCfg;1f;1000f);

greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

posResult:.commodity.modelreport.runPositionRisk[position;greekCfg;disCfg];

expectedKeys:`tradeId`commodity`comparison`disagreement`priceRows`greeksRows`scenarioPnlRows`alertRows`positionSummary`status`errorMessage;
.testutil.assertTrue[all expectedKeys in key posResult;"runPositionRisk has expected keys"];

.testutil.assertTrue[`T_PR=posResult`tradeId;"tradeId echoed"];
.testutil.assertTrue[`wti=posResult`commodity;"commodity echoed"];
.testutil.assertTrue[`OK=posResult`status;"position status OK"];

priceRowsTbl:posResult`priceRows;
.testutil.assertTrue[4=count priceRowsTbl;"priceRows has four model rows"];
.testutil.assertTrue[all (priceRowsTbl`tradeId)=`T_PR;"priceRows all tagged with tradeId"];
.testutil.assertTrue[all (priceRowsTbl`commodity)=`wti;"priceRows all tagged with commodity"];

greeksRowsTbl:posResult`greeksRows;
.testutil.assertTrue[4=count greeksRowsTbl;"greeksRows has four model rows"];
.testutil.assertTrue[all (greeksRowsTbl`tradeId)=`T_PR;"greeksRows all tagged with tradeId"];

alertRowsTbl:posResult`alertRows;
.testutil.assertTrue[6=count alertRowsTbl;"alertRows has six alert rows"];
.testutil.assertTrue[all (alertRowsTbl`tradeId)=`T_PR;"alertRows all tagged with tradeId"];

positionSummaryDict:posResult`positionSummary;
.testutil.assertTrue[`OK=positionSummaryDict`positionStatus;"positionSummary status OK"];
.testutil.assertTrue[1f=positionSummaryDict`quantity;"positionSummary quantity echoed"];
.testutil.assertTrue[1000f=positionSummaryDict`contractMultiplier;"positionSummary multiplier echoed"];

-1 "PASS test_modelreport_run_position_risk: priceRows=",string[count priceRowsTbl],", alertRows=",string[count alertRowsTbl];
