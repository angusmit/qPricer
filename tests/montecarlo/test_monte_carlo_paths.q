/ test_monte_carlo_paths.q - GBM path simulation
\l lib/init.q

mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(
    1000;50;42;0b;0b;0.95);
pathMatrix:.montecarlo.simulateGBMPaths[100f;0.05;0f;0.2;1f;mcConfig];

/ 1. Correct dimensions
.testutil.assertTrue[1000=count pathMatrix;"1000 paths"];
.testutil.assertTrue[50=count pathMatrix 0;"50 time steps per path"];

/ 2. All prices positive
allPositive:all {all x>0f} each pathMatrix;
.testutil.assertTrue[allPositive;"all simulated prices positive"];

/ 3. Same seed gives same paths
pathMatrix2:.montecarlo.simulateGBMPaths[100f;0.05;0f;0.2;1f;mcConfig];
.testutil.assertTrue[(pathMatrix 0)~pathMatrix2 0;"same seed gives same paths"];

/ 4. Different seed gives different paths
mcConfig2:@[mcConfig;`randomSeed;:;99];
pathMatrix3:.montecarlo.simulateGBMPaths[100f;0.05;0f;0.2;1f;mcConfig2];
.testutil.assertTrue[not (pathMatrix 0)~pathMatrix3 0;"different seed gives different paths"];

-1 "PASS test_monte_carlo_paths: paths=1000, steps=50";
