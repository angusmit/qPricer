\l lib/init.q
/ fromSynthetic produces a GBM path with the standard schema. Determinism with a
/ fixed seed must hold; different seeds must yield different paths; invalid configs
/ must be rejected; empty paths must fail validation.

pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0.05;0.30;20;1f%252f;0.05;0f;42);
pathTbl:.strategy.path.fromSynthetic pathCfg;

requiredCols:`stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status;
.testutil.assertTableColumns[pathTbl;requiredCols;"synthetic path schema"];
.testutil.assertTrue[20=count pathTbl;"length = steps"];
.testutil.assertTrue[100f=pathTbl[0;`spot];"first spot = spot0"];
.testutil.assertTrue[all (pathTbl`status)=`OK;"every status OK"];
.testutil.assertTrue[all 0f<pathTbl`spot;"every spot positive"];
.testutil.assertTrue[(til 20)~pathTbl`stepIndex;"stepIndex = 0..N-1"];

pathTbl2:.strategy.path.fromSynthetic pathCfg;
.testutil.assertTrue[(pathTbl`spot)~pathTbl2`spot;"deterministic with same seed"];

pathCfg3:@[pathCfg;`seed;:;99];
pathTbl3:.strategy.path.fromSynthetic pathCfg3;
.testutil.assertTrue[not (pathTbl`spot)~pathTbl3`spot;"different seed yields different path"];

zeroVolCfg:@[pathCfg;(`volatility;`drift);:;(0f;0f)];
flatPath:.strategy.path.fromSynthetic zeroVolCfg;
.testutil.assertTrue[all (flatPath`spot)=100f;"zero vol + zero drift yields flat path"];

badStepsResult:.[.strategy.path.fromSynthetic;enlist @[pathCfg;`steps;:;1];{`ERROR}];
.testutil.assertTrue[badStepsResult~`ERROR;"steps <= 1 rejected"];

emptyValResult:.[.strategy.path.validate;enlist ();{`ERROR}];
.testutil.assertTrue[emptyValResult~`ERROR;"empty pathTable rejected"];

-1 "PASS test_strategy_path_synthetic: rows=",string[count pathTbl],", lastSpot=",string last pathTbl`spot;
