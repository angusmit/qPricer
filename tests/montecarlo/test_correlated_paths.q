/ test_correlated_paths.q - correlated normal generation and multi-asset GBM
\l core/init.q

corrMatrix:(1 0.5 0.3;0.5 1 0.4;0.3 0.4 1f);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(
    5000;10;42;0b;0b;0.95);

/ 1. Correlated normals: correct dimensions
corrNormals:.montecarlo.generateCorrelatedNormals[5000;10;corrMatrix;42];
.testutil.assertTrue[50000=count corrNormals;"50000 rows (5000*10)"];
.testutil.assertTrue[3=count corrNormals 0;"3 symbols per row"];

/ 2. Correlated GBM paths
spotVector:100 250 800f;
rateVector:0.05 0.05 0.05;
divVector:0 0.01 0f;
volVector:0.2 0.25 0.35;
pathData:.montecarlo.simulateCorrelatedGBMPaths[spotVector;rateVector;divVector;volVector;1f;corrMatrix;mcConfig];
.testutil.assertTrue[3=count pathData;"3 symbol path sets"];
.testutil.assertTrue[5000=count pathData 0;"5000 paths for symbol 0"];
.testutil.assertTrue[10=count (pathData 0) 0;"10 steps per path"];

/ 3. All prices positive
allPositive:all {all {all x>0f} each x} each pathData;
.testutil.assertTrue[allPositive;"all simulated prices positive"];

/ 4. Same seed gives same paths
pathData2:.montecarlo.simulateCorrelatedGBMPaths[spotVector;rateVector;divVector;volVector;1f;corrMatrix;mcConfig];
.testutil.assertTrue[((pathData 0) 0)~(pathData2 0) 0;"same seed gives same paths"];

/ 5. Different seed gives different paths
mcConfig2:@[mcConfig;`randomSeed;:;99];
pathData3:.montecarlo.simulateCorrelatedGBMPaths[spotVector;rateVector;divVector;volVector;1f;corrMatrix;mcConfig2];
.testutil.assertTrue[not ((pathData 0) 0)~(pathData3 0) 0;"different seed gives different paths"];

-1 "PASS test_correlated_paths";
