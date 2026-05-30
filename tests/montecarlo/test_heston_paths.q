/ test_heston_paths.q
\l core/init.q
hestonParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;0.05;0.0);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;50;42;0b;0b;0.95);

pathResult:.heston.simulateHestonPaths[100f;hestonParams;1f;mcConfig];
.testutil.assertTrue[pathResult[`status]~`OK;"status OK"];
.testutil.assertTrue[5000=count pathResult`spotPaths;"5000 spot paths"];
.testutil.assertTrue[51=count (pathResult`spotPaths) 0;"51 steps (incl initial)"];
.testutil.assertTrue[5000=count pathResult`finalSpot;"5000 final spots"];
.testutil.assertTrue[all pathResult[`finalSpot]>0f;"all final spots positive"];

/ Same seed deterministic
pathResult2:.heston.simulateHestonPaths[100f;hestonParams;1f;mcConfig];
.testutil.assertTrue[(pathResult`finalSpot)~pathResult2`finalSpot;"same seed deterministic"];

/ Different seed different
mcConfig2:@[mcConfig;`randomSeed;:;99];
pathResult3:.heston.simulateHestonPaths[100f;hestonParams;1f;mcConfig2];
.testutil.assertTrue[not (pathResult`finalSpot)~pathResult3`finalSpot;"different seed different"];

/ Diagnostics
diagResult:.heston.pathDiagnostics pathResult;
.testutil.assertTrue[diagResult[`pathCount]=5000;"diag pathCount"];
.testutil.assertTrue[diagResult[`averageFinalSpot]>0f;"diag avgFinalSpot positive"];

-1 "PASS test_heston_paths: avgFinalSpot=",string[diagResult`averageFinalSpot],", avgFinalVar=",string[diagResult`averageFinalVariance],", negVarCount=",string diagResult`negativeRawVarianceCount;
