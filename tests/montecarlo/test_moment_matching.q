/ test_moment_matching.q
\l core/init.q
normalMatrix:.montecarlo.generateNormalPaths[5000;50;42];
mmMatrix:.variance.applyMomentMatching normalMatrix;

allNormals:raze mmMatrix;
sampleMean:avg allNormals;
sampleStd:dev allNormals;

.testutil.assertNear[sampleMean;0f;0.001;"moment matched mean ~ 0"];
.testutil.assertNear[sampleStd;1f;0.001;"moment matched std ~ 1"];

/ MC price with moment matching should still be close to BS
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(10000;1;42;0b;0b;0.95);
mcResult:.montecarlo.priceEuropeanMC[`call;100f;100f;1f;0.05;0f;0.2;mcConfig];
.testutil.assertTrue[(abs mcResult[`price]-bsCall)<1f;"MC with moment matching close to BS"];

-1 "PASS test_moment_matching: mean=",string[sampleMean],", std=",string sampleStd;
