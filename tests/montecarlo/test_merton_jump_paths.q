/ test_merton_jump_paths.q
\l lib/init.q
mertonParams:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.2;0.5;-0.1;0.3;0.05;0.0);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;50;42;0b;0b;0.95);

pathResult:.merton.simulateJumpDiffusionPaths[100f;mertonParams;1f;mcConfig];
.testutil.assertTrue[pathResult[`status]~`OK;"status OK"];
.testutil.assertTrue[5000=count pathResult`finalSpot;"5000 final spots"];
.testutil.assertTrue[all pathResult[`finalSpot]>0f;"all final spots positive"];
avgJumps:avg `float$pathResult`totalJumpCount;
.testutil.assertTrue[avgJumps>0f;"some jumps occurred"];

/ Deterministic
pathResult2:.merton.simulateJumpDiffusionPaths[100f;mertonParams;1f;mcConfig];
.testutil.assertTrue[(pathResult`finalSpot)~pathResult2`finalSpot;"same seed deterministic"];

/ Different seed
mcConfig2:@[mcConfig;`randomSeed;:;99];
pathResult3:.merton.simulateJumpDiffusionPaths[100f;mertonParams;1f;mcConfig2];
.testutil.assertTrue[not (pathResult`finalSpot)~pathResult3`finalSpot;"different seed different"];

diagResult:.merton.pathDiagnostics pathResult;
.testutil.assertTrue[diagResult[`pathCount]=5000;"diag pathCount"];

-1 "PASS test_merton_jump_paths: avgJumps=",string[avgJumps],", avgFinalSpot=",string diagResult`averageFinalSpot;
