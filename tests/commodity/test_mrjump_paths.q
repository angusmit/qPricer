\l lib/init.q
/ simulatePaths returns dict with the right shape and positive prices.

params:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0f;0.25);
x0:log 70f;
expiryVal:1f;

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;2000];
mcConfig:@[mcConfig;`timeStepCount;:;25];

pathResult:.commodity.mrjump.simulatePaths[x0;params;expiryVal;mcConfig];

.testutil.assertTrue[`pathMatrix in key pathResult;"result has pathMatrix"];
.testutil.assertTrue[`jumpCountMatrix in key pathResult;"result has jumpCountMatrix"];
.testutil.assertTrue[`jumpSizeMatrix in key pathResult;"result has jumpSizeMatrix"];

pathMatrix:pathResult`pathMatrix;
.testutil.assertTrue[2000=count pathMatrix;"pathCount rows"];
.testutil.assertTrue[25=count pathMatrix 0;"stepCount columns"];
.testutil.assertTrue[all all each pathMatrix>0f;"all spot prices positive"];

jumpCountMatrix:pathResult`jumpCountMatrix;
.testutil.assertTrue[2000=count jumpCountMatrix;"jumpCountMatrix has pathCount rows"];
.testutil.assertTrue[25=count jumpCountMatrix 0;"jumpCountMatrix has stepCount cols"];
.testutil.assertTrue[all all each jumpCountMatrix>=0;"all jump counts non-negative"];

pathDiag:.commodity.mrjump.pathDiagnostics pathMatrix;
.testutil.assertTrue[`pathCount in key pathDiag;"pathDiagnostics has pathCount"];
.testutil.assertTrue[`averageFinal in key pathDiag;"pathDiagnostics has averageFinal"];

-1 "PASS test_mrjump_paths: avgFinal=",string pathDiag`averageFinal;
