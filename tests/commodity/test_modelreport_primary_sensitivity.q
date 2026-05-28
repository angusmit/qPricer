\l lib/init.q
/ primarySensitivity maps each model row in greeksAllModels to its
/ per-model primary sensitivity:
/   black76   -> forwardDelta
/   schwartz  -> logStateDelta
/   schwartz2 -> longFactorSensitivity
/   mrjump    -> logStateDelta

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
primaryTbl:.commodity.modelreport.primarySensitivity greeksTbl;

requiredCols:`modelName`primarySensitivityName`primarySensitivity`pricingStatus`pricingErrorMessage;
.testutil.assertTableColumns[primaryTbl;requiredCols;"primarySensitivity schema"];
.testutil.assertTrue[4=count primaryTbl;"primarySensitivity has four rows"];

mappingByModel:exec primarySensitivityName by modelName from primaryTbl;
.testutil.assertTrue[`forwardDelta=first mappingByModel`black76;"black76 -> forwardDelta"];
.testutil.assertTrue[`logStateDelta=first mappingByModel`schwartz;"schwartz -> logStateDelta"];
.testutil.assertTrue[`longFactorSensitivity=first mappingByModel`schwartz2;"schwartz2 -> longFactorSensitivity"];
.testutil.assertTrue[`logStateDelta=first mappingByModel`mrjump;"mrjump -> logStateDelta"];

okMask:(primaryTbl`pricingStatus)=`OK;
okRows:primaryTbl where okMask;
.testutil.assertTrue[all not null okRows`primarySensitivity;"OK rows have non-null primarySensitivity"];

-1 "PASS test_modelreport_primary_sensitivity: rows=",string[count primaryTbl],", okCount=",string[count okRows];
