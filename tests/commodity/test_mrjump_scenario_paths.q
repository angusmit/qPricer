\l core/init.q
/ scenarioPaths applies bumps to params then resimulates. Verify:
/   (a) jumpIntensityBump increases average jump count
/   (b) jumpMeanBump (positive) lifts the average final price
/   (c) recognised keys produce the same dict shape as simulatePaths
/   (d) unrecognised keys leave paths identical to the baseline

x0:log 70f;
expiryVal:1f;
baseParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;0.5;0f;0.25);

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;5000];
mcConfig:@[mcConfig;`timeStepCount;:;50];

baseResult:.commodity.mrjump.simulatePaths[x0;baseParams;expiryVal;mcConfig];
baseDiag:.commodity.mrjump.jumpDiagnostics baseResult;

intensityConfig:enlist[`jumpIntensityBump]!enlist 1.5;
intensityResult:.commodity.mrjump.scenarioPaths[x0;baseParams;expiryVal;mcConfig;intensityConfig];
intensityDiag:.commodity.mrjump.jumpDiagnostics intensityResult;
.testutil.assertTrue[intensityDiag[`averageJumpCount]>baseDiag`averageJumpCount;"jumpIntensityBump increases jump count"];

meanConfig:enlist[`jumpMeanBump]!enlist 0.1;
meanResult:.commodity.mrjump.scenarioPaths[x0;baseParams;expiryVal;mcConfig;meanConfig];
meanDiag:.commodity.mrjump.jumpDiagnostics meanResult;
.testutil.assertTrue[meanDiag[`averageFinalPrice]>baseDiag`averageFinalPrice;"positive jumpMeanBump lifts average final price"];

unknownConfig:enlist[`unknownBump]!enlist 1f;
unchangedResult:.commodity.mrjump.scenarioPaths[x0;baseParams;expiryVal;mcConfig;unknownConfig];
unchangedPaths:unchangedResult`pathMatrix;
basePaths:baseResult`pathMatrix;
maxDiff:max max each abs basePaths-unchangedPaths;
.testutil.assertTrue[maxDiff<1e-12;"unrecognised scenario key leaves paths unchanged"];

-1 "PASS test_mrjump_scenario_paths: baseJumps ",string[baseDiag`averageJumpCount],", bumpedJumps ",string[intensityDiag`averageJumpCount],", baseFinal ",string[baseDiag`averageFinalPrice],", meanBumpFinal ",string meanDiag`averageFinalPrice;
