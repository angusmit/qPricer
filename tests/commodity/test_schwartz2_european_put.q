\l core/init.q
/ Closed-form European put vs Monte Carlo, plus exact put-call parity.

params:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.35;0.15;0f;0.3);
shortFactor0:0f;
longFactor0:log 75f;
strikeVal:75f;
expiryVal:1f;
rateVal:0.05;

closedFormCall:.commodity.schwartz2.europeanOptionPrice[`call;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;params];
closedFormPut:.commodity.schwartz2.europeanOptionPrice[`put;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;params];
.testutil.assertTrue[closedFormPut>0f;"put positive"];

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;1];
mcResult:.commodity.schwartz2.europeanOptionPriceMC[`put;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;params;mcConfig];
mcPrice:mcResult`price;
seVal:mcResult`standardError;
.testutil.assertTrue[(abs closedFormPut-mcPrice)<4f*seVal;"closed form put within 4 SE of MC"];

fwdValue:.commodity.schwartz2.futuresPrice[shortFactor0;longFactor0;params;expiryVal];
parityLhs:closedFormCall-closedFormPut;
parityRhs:(exp neg rateVal*expiryVal)*fwdValue-strikeVal;
.testutil.assertNear[parityLhs;parityRhs;1e-6;"put-call parity exp(-rT)*(F-K)"];

-1 "PASS test_schwartz2_european_put: put ",string[closedFormPut],", parity diff ",string abs parityLhs-parityRhs;
