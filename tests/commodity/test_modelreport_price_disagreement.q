\l lib/init.q
/ priceDisagreement summarises priceAllModels into min/max/range/pct and
/ raises priceRangeAlert when either the absolute or relative range exceeds
/ its threshold. Verify range arithmetic plus a threshold-trigger case.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;10000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

modelInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

priceTbl:.commodity.modelreport.priceAllModels[optionSetup;modelInputs];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];
priceDis:.commodity.modelreport.priceDisagreement[priceTbl;disCfg];

expectedKeys:`okModelCount`minPrice`maxPrice`averagePrice`priceRange`priceRangePct`priceRangeAlert`priceRangeAbsThreshold`priceRangePctThreshold`status`errorMessage;
.testutil.assertTrue[all expectedKeys in key priceDis;"priceDisagreement dict has expected keys"];
.testutil.assertTrue[4=priceDis`okModelCount;"okModelCount = 4"];
.testutil.assertTrue[`OK=priceDis`status;"status OK when all OK"];
.testutil.assertTrue[priceDis[`maxPrice]>=priceDis`minPrice;"max >= min"];
.testutil.assertNear[priceDis`priceRange;priceDis[`maxPrice]-priceDis`minPrice;1e-12;"priceRange = max - min"];
.testutil.assertNear[priceDis`priceRangePct;priceDis[`priceRange]%abs priceDis`averagePrice;1e-12;"priceRangePct = range / |average|"];

absThr:disCfg`priceRangeAbsThreshold;
pctThr:disCfg`priceRangePctThreshold;
absAlertExpected:priceDis[`priceRange]>absThr;
pctAlertExpected:(not null priceDis`priceRangePct)&priceDis[`priceRangePct]>pctThr;
expectedAlert:absAlertExpected|pctAlertExpected;
.testutil.assertTrue[priceDis[`priceRangeAlert]=expectedAlert;"alert flag matches threshold rule"];
.testutil.assertTrue[-1h=type priceDis`priceRangeAlert;"alert is boolean atom"];

strictCfg:@[disCfg;`priceRangeAbsThreshold;:;0.001];
strictCfg:@[strictCfg;`priceRangePctThreshold;:;0.001];
strictDis:.commodity.modelreport.priceDisagreement[priceTbl;strictCfg];
.testutil.assertTrue[strictDis`priceRangeAlert;"strict thresholds force alert"];

-1 "PASS test_modelreport_price_disagreement: range=",string[priceDis`priceRange],", pct=",string[priceDis`priceRangePct],", alert=",string[priceDis`priceRangeAlert];
