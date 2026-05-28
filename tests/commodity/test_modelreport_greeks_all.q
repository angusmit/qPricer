\l lib/init.q
/ greeksAllModels combines per-model greeks tables via union join, producing
/ a single table whose column set is the union of all model-specific schemas.
/ Cells for non-applicable sensitivities are null (e.g. black76 has no
/ jumpMeanSensitivity).

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
greeksTable:.commodity.modelreport.greeksAllModels[optionSetup;modelInputs;greekCfg];

.testutil.assertTrue[4=count greeksTable;"greeksAllModels has four rows"];

expectedModels:`black76`schwartz`schwartz2`mrjump;
actualModels:greeksTable`modelName;
.testutil.assertTrue[all expectedModels in actualModels;"all expected models present"];
.testutil.assertTrue[all actualModels in expectedModels;"no unexpected models"];

.testutil.assertTrue[all (greeksTable`pricingStatus)=`OK;"all model greeks status OK"];
.testutil.assertTrue[all (greeksTable`basePrice)>0f;"all basePrices positive"];

mergedCols:cols greeksTable;
.testutil.assertTrue[`forwardDelta in mergedCols;"merged table has black76 forwardDelta"];
.testutil.assertTrue[`logStateDelta in mergedCols;"merged table has logStateDelta from schwartz/mrjump"];
.testutil.assertTrue[`longFactorSensitivity in mergedCols;"merged table has schwartz2 longFactorSensitivity"];
.testutil.assertTrue[`jumpIntensitySensitivity in mergedCols;"merged table has mrjump jumpIntensitySensitivity"];
.testutil.assertTrue[`correlationSensitivity in mergedCols;"merged table has schwartz2 correlationSensitivity"];

black76RowSelect:greeksTable first where (greeksTable`modelName)=`black76;
.testutil.assertTrue[null black76RowSelect`jumpIntensitySensitivity;"black76 row has null jumpIntensitySensitivity"];
schwartz2RowSelect:greeksTable first where (greeksTable`modelName)=`schwartz2;
.testutil.assertTrue[null schwartz2RowSelect`forwardDelta;"schwartz2 row has null forwardDelta"];

-1 "PASS test_modelreport_greeks_all: rows=",string[count greeksTable],", cols=",string[count mergedCols];
