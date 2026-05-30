\l core/init.q
/ priceAllModels produces one row per registered commodity model.
/ All four prices are positive, all status flags OK, and the schema is fixed.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;20000];
mcCfg:@[mcCfg;`timeStepCount;:;50];

modelInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

priceTable:.commodity.modelreport.priceAllModels[optionSetup;modelInputs];

requiredCols:`modelName`modelPrice`standardError`lowerConfidence`upperConfidence`pricingStatus`pricingErrorMessage;
.testutil.assertTableColumns[priceTable;requiredCols;"priceAllModels has full schema"];
.testutil.assertTrue[4=count priceTable;"priceAllModels has four rows"];

expectedModels:`black76`schwartz`schwartz2`mrjump;
actualModels:priceTable`modelName;
.testutil.assertTrue[all expectedModels in actualModels;"expected models all present"];
.testutil.assertTrue[all actualModels in expectedModels;"no unexpected models"];

.testutil.assertTrue[all (priceTable`modelPrice)>0f;"all model prices positive"];
.testutil.assertTrue[all (priceTable`pricingStatus)=`OK;"all status OK"];

mrjumpRow:priceTable first where (priceTable`modelName)=`mrjump;
.testutil.assertTrue[mrjumpRow[`standardError]>0f;"mrjump standardError positive"];
.testutil.assertTrue[mrjumpRow[`upperConfidence]>mrjumpRow`lowerConfidence;"mrjump confidence band ordered"];

-1 "PASS test_modelreport_price_all: prices=",-3!priceTable`modelPrice;
