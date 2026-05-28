\l lib/init.q
/ jumpDiagnostics produces sensible jump statistics:
/   averageJumpCount ~ lambda * T
/   jumpPathFraction ~ 1 - exp(-lambda * T)

lambdaVal:2f;
expiryVal:1f;
params:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;lambdaVal;0.0;0.25);
x0:log 70f;

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;5000];
mcConfig:@[mcConfig;`timeStepCount;:;50];

pathResult:.commodity.mrjump.simulatePaths[x0;params;expiryVal;mcConfig];
diag:.commodity.mrjump.jumpDiagnostics pathResult;

expectedKeys:`pathCount`timeStepCount`averageFinalPrice`averageMinimumPrice`averageMaximumPrice`averageJumpCount`maxJumpCount`jumpPathFraction`status`errorMessage;
.testutil.assertTrue[all expectedKeys in key diag;"diagnostics dict has all expected keys"];
.testutil.assertTrue[diag[`status]=`OK;"diagnostics status OK"];

expectedAvgJumps:lambdaVal*expiryVal;
.testutil.assertNear[diag`averageJumpCount;expectedAvgJumps;0.1;"average jump count near lambda*T"];

expectedJumpFraction:1f-exp neg lambdaVal*expiryVal;
.testutil.assertNear[diag`jumpPathFraction;expectedJumpFraction;0.02;"jump path fraction near 1-exp(-lambda T)"];

.testutil.assertTrue[diag[`maxJumpCount]>=diag`averageJumpCount;"max jump count >= average"];
.testutil.assertTrue[diag[`averageMinimumPrice]<=diag`averageFinalPrice;"avg min <= avg final"];
.testutil.assertTrue[diag[`averageMaximumPrice]>=diag`averageFinalPrice;"avg max >= avg final"];

-1 "PASS test_mrjump_jump_diagnostics: avgJumps=",string[diag`averageJumpCount]," vs ",string[expectedAvgJumps],"; jumpFrac=",string[diag`jumpPathFraction]," vs ",string expectedJumpFraction;
