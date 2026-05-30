\l core/init.q
/ Verify simulated mean and variance of X_tau match exact OU formulas.

params:`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30);
x0:log 70f;
tauVal:0.5;

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;1];

pathMatrix:.commodity.schwartz.simulatePaths[x0;params;tauVal;mcConfig];
terminalLogs:log last each pathMatrix;
sampleMean:avg terminalLogs;
sampleVariance:var terminalLogs;

analytical:.commodity.schwartz.transitionMoments[x0;params;tauVal];
.testutil.assertNear[sampleMean;analytical`mean;0.01;"OU conditional mean matches"];
.testutil.assertNear[sampleVariance;analytical`variance;0.001;"OU conditional variance matches"];

-1 "PASS test_schwartz_transition_moments: mean ",string[sampleMean]," vs ",string[analytical`mean],"; var ",string[sampleVariance]," vs ",string analytical`variance;
