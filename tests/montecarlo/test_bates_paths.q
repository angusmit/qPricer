\l lib/init.q
batesParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;0.5;-0.1;0.3;0.05;0.0);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;50;42;0b;0b;0.95);
pathResult:.bates.simulateBatesPaths[100f;batesParams;1f;mcConfig];
.testutil.assertTrue[pathResult[`status]~`OK;"status OK"];
.testutil.assertTrue[5000=count pathResult`finalSpot;"5000 final spots"];
.testutil.assertTrue[all pathResult[`finalSpot]>0f;"all final spots positive"];
avgJumps:avg `float$pathResult`totalJumpCount;
.testutil.assertTrue[avgJumps>0f;"some jumps occurred"];
pathResult2:.bates.simulateBatesPaths[100f;batesParams;1f;mcConfig];
.testutil.assertTrue[(pathResult`finalSpot)~pathResult2`finalSpot;"same seed deterministic"];
diagResult:.bates.pathDiagnostics pathResult;
.testutil.assertTrue[diagResult[`pathCount]=5000;"diag pathCount"];
-1 "PASS test_bates_paths: avgFinalSpot=",string[diagResult`averageFinalSpot],", avgJumps=",string avgJumps;
