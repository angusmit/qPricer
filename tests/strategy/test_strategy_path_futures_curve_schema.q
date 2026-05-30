\l core/init.q
tenors:0.1 0.25 0.5 1f;
nSteps:5;
pathCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    tenors;
    `simple;
    `spot0`drift`volatility`contango!(50f;0.05;0.30;1.5);
    nSteps;
    1f%252f;
    21);
bundle:.strategy.path.fromFuturesCurve pathCfg;
requiredKeys:`frontPath`curveSnapshots`tenors`evolutionModel`frontLevels`jumpCountsAtStep;
.testutil.assertTrue[all requiredKeys in key bundle;"bundle keys present"];
frontPath:bundle`frontPath;
frontPathCols:`stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status;
.testutil.assertTableColumns[frontPath;frontPathCols;"frontPath standard schema"];
.testutil.assertTrue[nSteps=count frontPath;"frontPath rows = steps"];
curveSnapshots:bundle`curveSnapshots;
csCols:`stepIndex`stepDate`tenor`futuresPrice;
.testutil.assertTableColumns[curveSnapshots;csCols;"curveSnapshots schema"];
expectedRows:nSteps*count tenors;
.testutil.assertTrue[expectedRows=count curveSnapshots;"curveSnapshots rows = steps*tenors"];
/ Contango identity at step 0: F_T = spot0 + contango*T
spot0:50f;
contango:1.5;
step0Rows:curveSnapshots where (curveSnapshots`stepIndex)=0;
expectedPrices:spot0+contango*step0Rows`tenor;
.testutil.assertTrue[1e-10>max abs (step0Rows`futuresPrice)-expectedPrices;"contango identity at step 0"];
/ Front level at step 0 = spot0
.testutil.assertTrue[1e-10>abs (first bundle`frontLevels)-spot0;"front level at step 0 = spot0"];
-1 "PASS test_strategy_path_futures_curve_schema: contango identity holds";
