\l lib/init.q
/ Closed-form European call vs Monte Carlo.

params:`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30);
x0:log 70f;
strikeVal:72f;
expiryVal:0.5;
rateVal:0.05;

closedFormCall:.commodity.schwartz.europeanOptionPrice[`call;x0;strikeVal;expiryVal;rateVal;params];
.testutil.assertTrue[closedFormCall>0f;"call positive"];

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;1];
mcResult:.commodity.schwartz.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;params;mcConfig];
mcPrice:mcResult`price;
seVal:mcResult`standardError;

.testutil.assertTrue[(abs closedFormCall-mcPrice)<4f*seVal;"closed form within 4 SE of MC"];

-1 "PASS test_schwartz_european_call: closed ",string[closedFormCall]," vs MC ",string[mcPrice]," (SE ",string[seVal],")";
