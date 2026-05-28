\l lib/init.q
/ With jumpMean = 0 and jumpVolatility = 0, jumps have no price effect even at
/ positive jumpIntensity. Paths must equal the jumpIntensity = 0 case exactly,
/ because the diffusion RNG stream is unchanged by the jump consumption.

x0:log 70f;
baseParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;0f;0f;0f);
zeroSizeParams:@[baseParams;`jumpIntensity;:;2f];
expiryVal:1f;

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;1000];
mcConfig:@[mcConfig;`timeStepCount;:;25];

baseResult:.commodity.mrjump.simulatePaths[x0;baseParams;expiryVal;mcConfig];
zeroSizeResult:.commodity.mrjump.simulatePaths[x0;zeroSizeParams;expiryVal;mcConfig];

basePaths:baseResult`pathMatrix;
zeroSizePaths:zeroSizeResult`pathMatrix;

maxDiff:max max each abs basePaths-zeroSizePaths;
.testutil.assertTrue[maxDiff<1e-12;"zero jump size yields identical paths to zero intensity (max diff ",string[maxDiff],")"];

baseDiag:.commodity.mrjump.jumpDiagnostics baseResult;
zeroSizeDiag:.commodity.mrjump.jumpDiagnostics zeroSizeResult;
.testutil.assertTrue[baseDiag[`averageJumpCount]=0f;"base case has zero average jumps"];
.testutil.assertTrue[zeroSizeDiag[`averageJumpCount]>0f;"zero-size case still records positive jump counts"];
.testutil.assertNear[baseDiag`averageFinalPrice;zeroSizeDiag`averageFinalPrice;1e-12;"average final price unchanged"];

-1 "PASS test_mrjump_zero_jump_size: maxDiff ",string[maxDiff],", baseAvgJumps ",string[baseDiag`averageJumpCount],", zeroSizeAvgJumps ",string zeroSizeDiag`averageJumpCount;
