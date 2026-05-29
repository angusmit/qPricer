\l lib/init.q
/ Multi-commodity correlated curve adapter: two commodities (power, gas) with a
/ shared correlation, per-name contango, three tenors, 6 steps.
pathCfg:`names`correlationMatrix`spot0s`drifts`vols`contangos`tenors`steps`stepYears`riskFreeRate`seed!(
    `power`gas;
    (1 0.4f;0.4 1f);
    50 3.5f;
    0 0f;
    0.45 0.35f;
    0.5 0.1f;
    0.08 0.25 0.5f;
    6;
    1f%252f;
    0.03;
    77);
bundle:.strategy.path.fromCorrelatedCurves pathCfg;

/ Top-level bundle keys.
.testutil.assertTrue[all `names`tenors`correlationMatrix`vols`riskFreeRate`steps`stepYears`curves in key bundle;"top-level bundle keys present"];
.testutil.assertTrue[(bundle`names)~`power`gas;"names preserved in order"];

/ Per-name bundles are key-compatible with fromFuturesCurve output.
powerBundle:(bundle`curves)`power;
gasBundle:(bundle`curves)`gas;
expectedKeys:`frontPath`curveSnapshots`tenors`evolutionModel`evolutionParams`frontLevels`jumpCountsAtStep;
.testutil.assertTrue[expectedKeys~key powerBundle;"per-name bundle key-compatible with fromFuturesCurve"];

/ frontPath obeys the standard strategy path schema and is run-ready.
.strategy.path.validate powerBundle`frontPath;
.testutil.assertTrue[6=count powerBundle`frontPath;"power frontPath has 6 rows"];
.testutil.assertTrue[6=count gasBundle`frontPath;"gas frontPath has 6 rows"];

/ Front levels start at the given spot0s.
.testutil.assertNear[first powerBundle`frontLevels;50f;1e-10;"power front0 = spot0"];
.testutil.assertNear[first gasBundle`frontLevels;3.5f;1e-10;"gas front0 = spot0"];

/ Curve snapshots: contango is additive per tenor. At step 0 the 0.5y power
/ tenor price = front0 + contango*0.5 = 50 + 0.5*0.5 = 50.25.
snaps:powerBundle`curveSnapshots;
.testutil.assertTrue[all `stepIndex`stepDate`tenor`futuresPrice in cols snaps;"curveSnapshots schema"];
.testutil.assertTrue[18=count snaps;"6 steps * 3 tenors = 18 snapshot rows"];
farPower0:first exec futuresPrice from snaps where stepIndex=0, tenor=0.5;
.testutil.assertNear[farPower0;50.25;1e-10;"power 0.5y tenor = front0 + contango*tenor"];

/ riskFreeRate carried onto each frontPath.
.testutil.assertTrue[all 0.03=(powerBundle`frontPath)`riskFreeRate;"riskFreeRate on frontPath"];

-1 "PASS test_strategy_path_correlated_curves_schema: powerFront0=",(string first powerBundle`frontLevels),", gasFront0=",string first gasBundle`frontLevels;
