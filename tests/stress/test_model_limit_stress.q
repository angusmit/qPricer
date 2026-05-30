\l core/init.q
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;25;42;0b;0b;0.95);
cfg:enlist[`mcConfig]!enlist mcConfig;

checkResults:.modelcheck.runCoreModelLimitChecks cfg;
summaryResult:.limitcheck.summary checkResults;
.testutil.assertTrue[summaryResult[`checkCount]>0;"stress checks ran"];
.testutil.assertTrue[summaryResult[`passedCount]>0;"stress checks passed"];

statusCol:checkResults`status;
okCount:sum statusCol=`OK;
.testutil.assertTrue[okCount>0;"some OK rows"];

-1 "PASS test_model_limit_stress: checks=",string[summaryResult`checkCount],", passed=",string summaryResult`passedCount;
