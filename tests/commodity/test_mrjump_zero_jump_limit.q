\l core/init.q
/ With jumpIntensity = 0, mrjump reduces to schwartz one-factor:
/   (a) MC path matrix matches schwartz simulatePaths element-wise (same seed and antithetic),
/   (b) MC call price matches schwartz closed-form within 4 SE.

x0:log 70f;
kappaVal:1.5;
thetaVal:log 75f;
sigmaVal:0.30;
expiryVal:1f;
strikeVal:75f;
rateVal:0.05;

mrjumpParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(kappaVal;thetaVal;sigmaVal;0f;0f;0f);
schwartzParams:`meanReversionSpeed`longRunLogMean`volatility!(kappaVal;thetaVal;sigmaVal);

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;20000];
mcConfig:@[mcConfig;`timeStepCount;:;50];

mrjumpResult:.commodity.mrjump.simulatePaths[x0;mrjumpParams;expiryVal;mcConfig];
mrjumpPaths:mrjumpResult`pathMatrix;
schwartzPaths:.commodity.schwartz.simulatePaths[x0;schwartzParams;expiryVal;mcConfig];

maxPathDiff:max max each abs mrjumpPaths-schwartzPaths;
.testutil.assertTrue[maxPathDiff<1e-12;"pathwise equality with schwartz (max diff = ",string[maxPathDiff],")"];

mcCallResult:.commodity.mrjump.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;mrjumpParams;mcConfig];
mcCallPrice:mcCallResult`price;
seVal:mcCallResult`standardError;
closedFormCall:.commodity.schwartz.europeanOptionPrice[`call;x0;strikeVal;expiryVal;rateVal;schwartzParams];

.testutil.assertTrue[(abs mcCallPrice-closedFormCall)<4f*seVal;"mrjump MC call within 4 SE of schwartz closed-form"];

-1 "PASS test_mrjump_zero_jump_limit: path diff ",string[maxPathDiff],", call MC ",string[mcCallPrice]," vs closed ",string[closedFormCall]," (SE ",string[seVal],")";
