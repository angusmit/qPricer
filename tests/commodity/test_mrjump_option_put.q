\l lib/init.q
/ MC European put is positive. Put-call parity at the sample level holds exactly:
/   C - P = exp(-rT) * (sampleMean(S_T) - K)
/ because both options price off the same MC sample (deterministic by seed).

x0:log 70f;
strikeVal:75f;
expiryVal:1f;
rateVal:0.05;
params:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0f;0.25);

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;20000];
mcConfig:@[mcConfig;`timeStepCount;:;50];

callResult:.commodity.mrjump.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;params;mcConfig];
putResult:.commodity.mrjump.europeanOptionPriceMC[`put;x0;strikeVal;expiryVal;rateVal;params;mcConfig];
.testutil.assertTrue[putResult[`price]>0f;"put price positive"];

pathResult:.commodity.mrjump.simulatePaths[x0;params;expiryVal;mcConfig];
sampleMean:avg last each pathResult`pathMatrix;
discFactor:exp neg rateVal*expiryVal;
parityLhs:(callResult`price)-putResult`price;
parityRhs:discFactor*sampleMean-strikeVal;
.testutil.assertNear[parityLhs;parityRhs;1e-10;"sample put-call parity"];

-1 "PASS test_mrjump_option_put: put ",string[putResult`price],", parity diff ",string abs parityLhs-parityRhs;
