\l core/init.q
/ scenarioShift applies forward / spot / vol bumps to the model input dictionary
/ and reprices all models; scenarioPnL produces per-model PnL relative to base.
/ Verify that scenarioPnL = quantity * contractMultiplier * (scenarioPrice - basePrice)
/ for every model row and that all expected models are represented.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;15000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

modelInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

scenarioConfig:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);

quantityVal:1f;
contractMultiplierVal:1000f;

basePriceTable:.commodity.modelreport.priceAllModels[optionSetup;modelInputs];
scenarioPriceTable:.commodity.modelreport.scenarioShift[optionSetup;modelInputs;scenarioConfig];
pnLTable:.commodity.modelreport.scenarioPnL[basePriceTable;scenarioPriceTable;quantityVal;contractMultiplierVal];

requiredCols:`modelName`basePrice`scenarioPrice`scenarioPnL`quantity`contractMultiplier`status`errorMessage;
.testutil.assertTableColumns[pnLTable;requiredCols;"scenarioPnL has full schema"];
.testutil.assertTrue[4=count pnLTable;"scenarioPnL has four rows"];
.testutil.assertTrue[all (pnLTable`status)=`OK;"all pnL rows OK"];
.testutil.assertTrue[all (pnLTable`quantity)=quantityVal;"quantity propagated"];
.testutil.assertTrue[all (pnLTable`contractMultiplier)=contractMultiplierVal;"contractMultiplier propagated"];

expectedPnL:contractMultiplierVal*quantityVal*((pnLTable`scenarioPrice)-pnLTable`basePrice);
maxPnLDiff:max abs (pnLTable`scenarioPnL)-expectedPnL;
.testutil.assertTrue[maxPnLDiff<1e-9;"scenarioPnL equals mult * qty * (scenario - base)"];

black76Row:pnLTable first where (pnLTable`modelName)=`black76;
black76Expected:contractMultiplierVal*quantityVal*((black76Row`scenarioPrice)-black76Row`basePrice);
.testutil.assertNear[black76Row`scenarioPnL;black76Expected;1e-9;"black76 scenarioPnL matches formula"];

-1 "PASS test_modelreport_scenario_pnl: maxFormulaDiff=",string[maxPnLDiff],", black76PnL=",string black76Row`scenarioPnL;
