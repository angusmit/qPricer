\l lib/init.q
/ Seasonality is OPT-IN: when seasonCfg is absent the curve adapters must be
/ byte-identical to the pre-v0.51 baseline (no extra keys, frontPath untouched).
base:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.04 0.25 0.5 1.0f;`simple;
    `spot0`drift`volatility`contango`riskFreeRate!(60f;0f;0.2;0.5;0.02);
    6;1f%252f;11);
b0:.strategy.path.fromFuturesCurve base;
b0Repeat:.strategy.path.fromFuturesCurve base;
.testutil.assertTrue[b0~b0Repeat;"fromFuturesCurve baseline deterministic / byte-identical"];

withSeason:base,enlist[`seasonCfg]!enlist `amplitude`phaseYears!(0.15;0f);
b1:.strategy.path.fromFuturesCurve withSeason;
.testutil.assertTrue[(key b0)~key b1;"seasonCfg adds NO bundle keys"];
.testutil.assertTrue[(b0`frontPath)~b1`frontPath;"frontPath unchanged by seasonality"];
.testutil.assertTrue[not (b0`curveSnapshots)~b1`curveSnapshots;"seasonality DOES alter curveSnapshots"];

/ fromCorrelatedCurves: same opt-in guarantee per name.
cc:`names`correlationMatrix`spot0s`drifts`vols`contangos`tenors`steps`stepYears`riskFreeRate`seed!(
    `a`b;(1 0.3f;0.3 1f);60 3f;0 0f;0.2 0.3f;0 0f;0.04 0.5 1.0f;6;1f%252f;0.02;9);
c0:.strategy.path.fromCorrelatedCurves cc;
.testutil.assertTrue[c0~.strategy.path.fromCorrelatedCurves cc;"fromCorrelatedCurves baseline deterministic"];
ccSeason:cc,enlist[`seasonCfg]!enlist `amplitude`phaseYears!(0.15;0f);
c1:.strategy.path.fromCorrelatedCurves ccSeason;
.testutil.assertTrue[(key c0)~key c1;"correlated seasonCfg adds NO bundle keys"];
.testutil.assertTrue[(c0[`curves;`a]`frontPath)~c1[`curves;`a]`frontPath;"correlated frontPath unchanged"];
.testutil.assertTrue[not (c0[`curves;`a]`curveSnapshots)~c1[`curves;`a]`curveSnapshots;"correlated seasonality alters curve a"];

-1 "PASS test_seasonality_overlay_baseline";
