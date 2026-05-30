\l core/init.q
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;25;42;0b;0b;0.95);
cfg:enlist[`mcConfig]!enlist mcConfig;

checkResults:.modelcheck.runCoreModelLimitChecks cfg;
.testutil.assertTrue[(count checkResults)>0;"checks returned rows"];

summaryResult:.limitcheck.summary checkResults;
.testutil.assertTrue[summaryResult[`checkCount]>0;"summary checkCount > 0"];
.testutil.assertTrue[summaryResult[`passedCount]>0;"at least one check passed"];

consistencyReport:.modelcheck.modelConsistencyReport checkResults;
.testutil.assertTrue[(count consistencyReport)>0;"consistency report has rows"];
.testutil.assertTrue[`severity in cols consistencyReport;"has severity column"];

-1 "PASS test_model_consistency_report: checks=",string[summaryResult`checkCount],", passed=",string[summaryResult`passedCount],", failed=",string summaryResult`failedCount;
