\l core/init.q
/ Positive jumpMean lifts the conditional mean of the terminal log price and
/ therefore the sample average final spot and a moderately OTM call price.
/ Compare jumpMean = 0 vs jumpMean = 0.10 with everything else fixed.

x0:log 70f;
strikeVal:75f;
expiryVal:1f;
rateVal:0.05;
zeroMeanParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;2f;0f;0.20);
positiveMeanParams:@[zeroMeanParams;`jumpMean;:;0.10];

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;50];

zeroMeanPathResult:.commodity.mrjump.simulatePaths[x0;zeroMeanParams;expiryVal;mcConfig];
positiveMeanPathResult:.commodity.mrjump.simulatePaths[x0;positiveMeanParams;expiryVal;mcConfig];

zeroMeanDiag:.commodity.mrjump.jumpDiagnostics zeroMeanPathResult;
positiveMeanDiag:.commodity.mrjump.jumpDiagnostics positiveMeanPathResult;

.testutil.assertNear[zeroMeanDiag`averageJumpCount;positiveMeanDiag`averageJumpCount;1e-12;"jump counts independent of jumpMean (same seed)"];
.testutil.assertTrue[positiveMeanDiag[`averageFinalPrice]>zeroMeanDiag`averageFinalPrice;"positive jumpMean lifts average final price"];

zeroMeanCallResult:.commodity.mrjump.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;zeroMeanParams;mcConfig];
positiveMeanCallResult:.commodity.mrjump.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;positiveMeanParams;mcConfig];
.testutil.assertTrue[positiveMeanCallResult[`price]>zeroMeanCallResult`price;"positive jumpMean lifts MC call price"];

-1 "PASS test_mrjump_jump_mean_sensitivity: zeroMeanFinal=",string[zeroMeanDiag`averageFinalPrice],", positiveMeanFinal=",string[positiveMeanDiag`averageFinalPrice],", zeroMeanCall=",string[zeroMeanCallResult`price],", positiveMeanCall=",string positiveMeanCallResult`price;
