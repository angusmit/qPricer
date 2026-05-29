\l lib/init.q
/ Determinism + correlation structure of the multi-commodity curve adapter.
baseCfg:`names`correlationMatrix`spot0s`drifts`vols`contangos`tenors`steps`stepYears`riskFreeRate`seed!(
    `crude`product;
    (1 0.7f;0.7 1f);
    80 95f;
    0 0f;
    0.30 0.35f;
    0.2 0.3f;
    0.08 0.25f;
    8;
    1f%252f;
    0.02;
    101);
bundleA:.strategy.path.fromCorrelatedCurves baseCfg;
bundleB:.strategy.path.fromCorrelatedCurves baseCfg;
/ Same seed -> byte-identical front levels for every name.
.testutil.assertTrue[(bundleA[`curves;`crude]`frontLevels)~bundleB[`curves;`crude]`frontLevels;"crude paths deterministic for fixed seed"];
.testutil.assertTrue[(bundleA[`curves;`product]`frontLevels)~bundleB[`curves;`product]`frontLevels;"product paths deterministic for fixed seed"];

/ Different seed -> different paths.
diffCfg:@[baseCfg;`seed;:;999];
bundleC:.strategy.path.fromCorrelatedCurves diffCfg;
.testutil.assertTrue[not (bundleA[`curves;`crude]`frontLevels)~bundleC[`curves;`crude]`frontLevels;"different seed -> different crude path"];

/ Correlation: the two commodities' log returns should be positively correlated
/ (input rho = 0.7). Recover realized correlation from the simulated fronts.
crudeRet:1_ deltas log bundleA[`curves;`crude]`frontLevels;
productRet:1_ deltas log bundleA[`curves;`product]`frontLevels;
realizedCorr:crudeRet cor productRet;
.testutil.assertTrue[realizedCorr>0.3;"realized correlation positive (rho=0.7 input)"];

/ Validation: non-symmetric correlation matrix rejected.
badCfg:@[baseCfg;`correlationMatrix;:;(1 0.7f;0.5 1f)];
badRes:@[.strategy.path.fromCorrelatedCurves;badCfg;{`ERROR}];
.testutil.assertTrue[badRes~`ERROR;"non-symmetric correlation matrix rejected"];

-1 "PASS test_strategy_path_correlated_curves_determinism: realizedCorr=",string realizedCorr;
