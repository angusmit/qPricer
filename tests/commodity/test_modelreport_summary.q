\l core/init.q
/ comparisonSummary returns modelCount / okModelCount / min / max / average /
/ price range plus status. Sanity-check ordering and counts under the standard
/ four-model setup (all-OK case).

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
summaryResult:.commodity.modelreport.comparisonSummary priceTable;

expectedKeys:`baselineModel`baselinePrice`averageModelPrice`modelCount`okModelCount`minPrice`maxPrice`priceRange`status`errorMessage;
.testutil.assertTrue[all expectedKeys in key summaryResult;"summary dict has all expected keys"];

.testutil.assertTrue[4=summaryResult`modelCount;"modelCount = 4"];
.testutil.assertTrue[4=summaryResult`okModelCount;"okModelCount = 4 (all priced)"];
.testutil.assertTrue[`OK=summaryResult`status;"status OK when all models priced"];

minPrice:summaryResult`minPrice;
avgPrice:summaryResult`averageModelPrice;
maxPrice:summaryResult`maxPrice;
priceRangeVal:summaryResult`priceRange;
.testutil.assertTrue[minPrice<=avgPrice;"min <= average"];
.testutil.assertTrue[avgPrice<=maxPrice;"average <= max"];
.testutil.assertTrue[priceRangeVal>=0f;"priceRange non-negative"];
.testutil.assertNear[priceRangeVal;maxPrice-minPrice;1e-12;"priceRange = max - min"];

.testutil.assertTrue[`black76=summaryResult`baselineModel;"baselineModel is black76"];
baselineRowSelect:priceTable first where (priceTable`modelName)=`black76;
.testutil.assertNear[summaryResult`baselinePrice;baselineRowSelect`modelPrice;1e-12;"baselinePrice matches priceTable black76 row"];

-1 "PASS test_modelreport_summary: count=",string[summaryResult`okModelCount]," min=",string[minPrice]," avg=",string[avgPrice]," max=",string maxPrice;
