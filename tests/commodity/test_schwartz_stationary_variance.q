\l lib/init.q
/ At long T the sample variance of X_T should approach sigma^2/(2 kappa).

params:`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30);
x0:log 70f;
expiryVal:10f;

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;50];

pathMatrix:.commodity.schwartz.simulatePaths[x0;params;expiryVal;mcConfig];
terminalLogs:log last each pathMatrix;
sampleVariance:var terminalLogs;
analyticalStationary:.commodity.schwartz.stationaryVariance params;

.testutil.assertNear[sampleVariance;analyticalStationary;0.002;"long-T variance approaches stationary"];

-1 "PASS test_schwartz_stationary_variance: sample ",string[sampleVariance]," vs stationary ",string analyticalStationary;
