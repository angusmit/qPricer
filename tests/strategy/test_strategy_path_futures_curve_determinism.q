\l lib/init.q
pathCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.1 0.25 0.5 1f;
    `simple;
    `spot0`drift`volatility`contango!(50f;0f;0.30;2f);
    6;
    1f%252f;
    99);
bundle1:.strategy.path.fromFuturesCurve pathCfg;
bundle2:.strategy.path.fromFuturesCurve pathCfg;
.testutil.assertTrue[(bundle1`frontLevels)~bundle2`frontLevels;"front levels deterministic"];
.testutil.assertTrue[(bundle1`curveSnapshots)~bundle2`curveSnapshots;"curve snapshots deterministic"];
.testutil.assertTrue[(bundle1`frontPath)~bundle2`frontPath;"front path deterministic"];
/ Different seed -> different paths
diffCfg:@[pathCfg;`seed;:;7];
bundle3:.strategy.path.fromFuturesCurve diffCfg;
.testutil.assertTrue[not (bundle3`frontLevels)~bundle1`frontLevels;"different seed -> different front"];
-1 "PASS test_strategy_path_futures_curve_determinism: simple model, 4 tenors x 6 steps";
