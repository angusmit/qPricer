\l core/init.q
/ Closed-form European put vs Monte Carlo and put-call parity.

params:`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30);
x0:log 70f;
strikeVal:72f;
expiryVal:0.5;
rateVal:0.05;

closedFormCall:.commodity.schwartz.europeanOptionPrice[`call;x0;strikeVal;expiryVal;rateVal;params];
closedFormPut:.commodity.schwartz.europeanOptionPrice[`put;x0;strikeVal;expiryVal;rateVal;params];
.testutil.assertTrue[closedFormPut>0f;"put positive"];

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;1];
mcResult:.commodity.schwartz.europeanOptionPriceMC[`put;x0;strikeVal;expiryVal;rateVal;params;mcConfig];
mcPrice:mcResult`price;
seVal:mcResult`standardError;
.testutil.assertTrue[(abs closedFormPut-mcPrice)<4f*seVal;"put closed form within 4 SE of MC"];

fwdValue:.commodity.schwartz.futuresPrice[x0;params;expiryVal];
parityLhs:closedFormCall-closedFormPut;
parityRhs:(exp neg rateVal*expiryVal)*fwdValue-strikeVal;
.testutil.assertNear[parityLhs;parityRhs;1e-6;"put-call parity exp(-rT)*(F-K)"];

-1 "PASS test_schwartz_european_put: put ",string[closedFormPut],", parity diff ",string abs parityLhs-parityRhs;
