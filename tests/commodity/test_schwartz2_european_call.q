\l lib/init.q
/ Closed-form European call (Black-76 on lognormal F, effective vol) vs Monte Carlo.

params:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.35;0.15;0f;0.3);
shortFactor0:0f;
longFactor0:log 75f;
strikeVal:75f;
expiryVal:1f;
rateVal:0.05;

closedFormCall:.commodity.schwartz2.europeanOptionPrice[`call;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;params];
.testutil.assertTrue[closedFormCall>0f;"call positive"];

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;1];
mcResult:.commodity.schwartz2.europeanOptionPriceMC[`call;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;params;mcConfig];
mcPrice:mcResult`price;
seVal:mcResult`standardError;

.testutil.assertTrue[(abs closedFormCall-mcPrice)<4f*seVal;"closed form call within 4 SE of MC"];

-1 "PASS test_schwartz2_european_call: closed ",string[closedFormCall]," vs MC ",string[mcPrice]," (SE ",string[seVal],")";
