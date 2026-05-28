/ test_mc_convergence_table.q
\l lib/init.q
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];

/ Create pricing projection
priceFn:{[pathCountVal]
    localCfg:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(pathCountVal;1;42;0b;0b;0.95);
    .montecarlo.priceEuropeanMC[`call;100f;100f;1f;0.05;0f;0.2;localCfg]
 };

convergenceTable:.convergence.mcConvergenceTable[priceFn;1000 5000 10000 25000;bsCall];

/ One row per pathCount
.testutil.assertTrue[4=count convergenceTable;"4 convergence rows"];

/ Error should generally decrease
errors:{x`absoluteError} each convergenceTable;
firstError:errors 0;
lastError:last errors;

/ SE should generally decrease
stdErrors:{x`standardError} each convergenceTable;
firstSE:stdErrors 0;
lastSE:last stdErrors;
.testutil.assertTrue[lastSE<=firstSE;"SE decreases with more paths"];

/ Summary
summaryResult:.convergence.summariseConvergence convergenceTable;
.testutil.assertTrue[summaryResult[`minPathCount]=1000;"min pathCount"];
.testutil.assertTrue[summaryResult[`maxPathCount]=25000;"max pathCount"];

-1 "PASS test_mc_convergence_table: firstError=",string[firstError],", lastError=",string[lastError],", firstSE=",string[firstSE],", lastSE=",string lastSE;
