\l lib/init.q
/ Verify simulated factor means, variances, and cross-covariance match analytics.

params:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.35;0.15;0f;0.3);
shortFactor0:0f;
longFactor0:log 75f;
tauVal:0.5;

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;1];

factorsResult:.commodity.schwartz2.simulateFactors[shortFactor0;longFactor0;params;tauVal;mcConfig];
shortPaths:factorsResult`shortPaths;
longPaths:factorsResult`longPaths;
shortTerminals:last each shortPaths;
longTerminals:last each longPaths;

sampleShortMean:avg shortTerminals;
sampleLongMean:avg longTerminals;
sampleShortVariance:var shortTerminals;
sampleLongVariance:var longTerminals;
sampleCrossCovariance:avg (shortTerminals-sampleShortMean)*longTerminals-sampleLongMean;

analytical:.commodity.schwartz2.factorMoments[shortFactor0;longFactor0;params;tauVal];

.testutil.assertNear[sampleShortMean;analytical`shortMean;0.01;"short factor mean matches"];
.testutil.assertNear[sampleLongMean;analytical`longMean;0.01;"long factor mean matches"];
.testutil.assertNear[sampleShortVariance;analytical`shortVariance;0.001;"short factor variance matches"];
.testutil.assertNear[sampleLongVariance;analytical`longVariance;0.001;"long factor variance matches"];
.testutil.assertNear[sampleCrossCovariance;analytical`crossCovariance;0.001;"cross covariance matches"];

-1 "PASS test_schwartz2_factor_moments: shortMean ",string[sampleShortMean]," vs ",string[analytical`shortMean],"; cov ",string[sampleCrossCovariance]," vs ",string analytical`crossCovariance;
