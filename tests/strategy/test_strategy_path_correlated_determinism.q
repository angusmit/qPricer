\l lib/init.q
pathCfg:`names`weights`spot0`vols`correlationMatrix`steps`stepYears`riskFreeRate`dividendYields`seed!(
    `A`B;0.5 0.5;100 80f;0.20 0.30;(1 0.4f;0.4 1f);6;1f%252f;0.02;0 0f;42);
bundle1:.strategy.path.fromCorrelated pathCfg;
bundle2:.strategy.path.fromCorrelated pathCfg;
namesA:bundle1`names;
.testutil.assertTrue[namesA~bundle2`names;"names deterministic"];
pathA1:(bundle1`pathTables)`A;
pathA2:(bundle2`pathTables)`A;
.testutil.assertTrue[(pathA1`spot)~pathA2`spot;"name A spot path deterministic"];
.testutil.assertTrue[((bundle1`pathTables)[`B]`spot)~(bundle2`pathTables)[`B]`spot;"name B spot path deterministic"];
.testutil.assertTrue[(bundle1`indexPath)~bundle2`indexPath;"indexPath deterministic"];
.testutil.assertTrue[(bundle1`correlationMatrix)~bundle2`correlationMatrix;"corr matrix returned"];
/ Different seed -> different paths
diffSeedCfg:@[pathCfg;`seed;:;99];
bundle3:.strategy.path.fromCorrelated diffSeedCfg;
.testutil.assertTrue[not (bundle3`indexPath)~bundle1`indexPath;"different seed -> different indexPath"];
-1 "PASS test_strategy_path_correlated_determinism: 2 names, 6 steps, seed=42 reproducible";
