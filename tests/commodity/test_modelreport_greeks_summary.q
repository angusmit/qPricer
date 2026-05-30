\l core/init.q
/ greeksSummary aggregates the merged greeks table into modelCount /
/ okModelCount / max-magnitude statistics, picking each model's primary
/ delta (forwardDelta for black76, logStateDelta for schwartz/mrjump,
/ longFactorSensitivity for schwartz2). When one model row is missing
/ inputs and fails, the summary must still report OK and reduce okModelCount.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;8000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

fullInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

greekCfg:.commodity.modelreport.defaultGreekConfig[];
greeksTable:.commodity.modelreport.greeksAllModels[optionSetup;fullInputs;greekCfg];
summaryDict:.commodity.modelreport.greeksSummary greeksTable;

expectedKeys:`modelCount`okModelCount`maxPrimaryDelta`maxVolatilityVega`maxJumpIntensitySensitivity`status`errorMessage;
.testutil.assertTrue[all expectedKeys in key summaryDict;"summary has expected keys"];
.testutil.assertTrue[4=summaryDict`modelCount;"modelCount = 4"];
.testutil.assertTrue[4=summaryDict`okModelCount;"okModelCount = 4 when all OK"];
.testutil.assertTrue[`OK=summaryDict`status;"status OK in all-OK case"];
.testutil.assertTrue[summaryDict[`maxPrimaryDelta]>0f;"maxPrimaryDelta positive"];
.testutil.assertTrue[summaryDict[`maxVolatilityVega]>0f;"maxVolatilityVega positive"];
.testutil.assertTrue[not null summaryDict`maxJumpIntensitySensitivity;"maxJumpIntensitySensitivity finite"];

partialInputs:`black76`schwartz!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30)));
partialTable:.commodity.modelreport.greeksAllModels[optionSetup;partialInputs;greekCfg];
partialSummary:.commodity.modelreport.greeksSummary partialTable;
.testutil.assertTrue[4=partialSummary`modelCount;"partial: modelCount counts all rows"];
.testutil.assertTrue[2=partialSummary`okModelCount;"partial: okModelCount counts OK rows only"];
.testutil.assertTrue[`OK=partialSummary`status;"partial: status OK with at least one OK model"];

emptyInputs:()!();
emptyTable:.commodity.modelreport.greeksAllModels[optionSetup;emptyInputs;greekCfg];
emptySummary:.commodity.modelreport.greeksSummary emptyTable;
.testutil.assertTrue[`ERROR=emptySummary`status;"empty: status ERROR when zero models OK"];
.testutil.assertTrue[0=emptySummary`okModelCount;"empty: okModelCount = 0"];

-1 "PASS test_modelreport_greeks_summary: maxPrimaryDelta=",string[summaryDict`maxPrimaryDelta]," maxVolVega=",string[summaryDict`maxVolatilityVega]," partialOk=",string[partialSummary`okModelCount]," / ",string[partialSummary`modelCount];
