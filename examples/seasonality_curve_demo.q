\l lib/init.q
/ ============================================================================
/ seasonality_curve_demo.q - synthetic demo of the opt-in seasonality overlay
/ (NOT a test, no real data). Crude is only mildly seasonal; this overlay is
/ intended for gas / power (e.g. natural gas, NG) curves once that data is added.
/ Here a gas-like front (~$3) is shaped by a winter-peaking annual cycle.
/ ============================================================================
baseCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.08 0.33 0.58 0.83 1.08f;`simple;
    `spot0`drift`volatility`contango`riskFreeRate!(3.0;0f;0.4;0.0;0.02);
    6;1f%252f;7);
flat:.strategy.path.fromFuturesCurve baseCfg;

winterPhase:`amplitude`phaseYears!(0.25;0f);
seasonal:.strategy.path.fromFuturesCurve baseCfg,enlist[`seasonCfg]!enlist winterPhase;

-1 "Baseline step-0 curve (no seasonality, flat at front level):";
show select tenor,futuresPrice from flat`curveSnapshots where stepIndex=0;
-1 "";
-1 "Seasonal step-0 curve (amplitude 0.25, winter-peaking cycle):";
show select tenor,futuresPrice from seasonal`curveSnapshots where stepIndex=0;
-1 "";
-1 "Note: tenors near integer years (winter) are shaped UP, half-year (summer) DOWN.";
