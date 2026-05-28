\l lib/init.q
/ priceDifferences reports per-model differenceToBaseline / pctDifferenceToBaseline
/ / differenceToAverage against a chosen baseline model. The baseline row itself
/ must have a zero difference and a zero percentage difference.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;15000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

modelInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

priceTable:.commodity.modelreport.priceAllModels[optionSetup;modelInputs];
diffsTable:.commodity.modelreport.priceDifferences[priceTable;`black76];

requiredCols:`modelName`modelPrice`baselineModel`baselinePrice`differenceToBaseline`pctDifferenceToBaseline`differenceToAverage`pricingStatus`pricingErrorMessage;
.testutil.assertTableColumns[diffsTable;requiredCols;"priceDifferences has full schema"];
.testutil.assertTrue[4=count diffsTable;"priceDifferences has four rows"];
.testutil.assertTrue[all (diffsTable`baselineModel)=`black76;"baselineModel populated everywhere"];

baselineRow:diffsTable first where (diffsTable`modelName)=`black76;
.testutil.assertNear[baselineRow`differenceToBaseline;0f;1e-12;"black76 differenceToBaseline = 0"];
.testutil.assertNear[baselineRow`pctDifferenceToBaseline;0f;1e-12;"black76 pctDifferenceToBaseline = 0"];

schwartzRow:diffsTable first where (diffsTable`modelName)=`schwartz;
.testutil.assertNear[schwartzRow`differenceToBaseline;schwartzRow[`modelPrice]-baselineRow`modelPrice;1e-12;"schwartz diff equals modelPrice minus baselinePrice"];

summary:.commodity.modelreport.comparisonSummary priceTable;
avgPrice:summary`averageModelPrice;
.testutil.assertNear[schwartzRow`differenceToAverage;schwartzRow[`modelPrice]-avgPrice;1e-12;"differenceToAverage equals modelPrice minus average"];

-1 "PASS test_modelreport_price_differences: schwartzDiff=",string[schwartzRow`differenceToBaseline],", schwartzPct=",string schwartzRow`pctDifferenceToBaseline;
