\l lib/init.q
/ modelRanking returns one row per OK model, sorted by descending modelPrice,
/ with rank 1 .. N. Verify ordering and that the row count matches the number
/ of OK rows in the input price table.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;10000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

modelInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

priceTable:.commodity.modelreport.priceAllModels[optionSetup;modelInputs];
okRowCount:count priceTable where (priceTable`pricingStatus)=`OK;
rankTable:.commodity.modelreport.modelRanking priceTable;

requiredCols:`rank`modelName`modelPrice;
.testutil.assertTableColumns[rankTable;requiredCols;"modelRanking has rank/modelName/modelPrice"];
.testutil.assertTrue[okRowCount=count rankTable;"rankTable row count matches okRowCount"];
.testutil.assertTrue[(rankTable`rank)~1+til count rankTable;"ranks are 1..N in order"];

modelPriceSeries:rankTable`modelPrice;
descendingPairs:(-1_modelPriceSeries)>=1_modelPriceSeries;
.testutil.assertTrue[all descendingPairs;"modelPrice non-increasing across ranks"];

-1 "PASS test_modelreport_ranking: count=",string[count rankTable],", top=",string[first modelPriceSeries],", bottom=",string last modelPriceSeries;
