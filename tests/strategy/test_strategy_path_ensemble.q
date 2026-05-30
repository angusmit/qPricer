\l core/init.q
/ path.ensemble produces deterministic numPaths paths via baseSeed+i seeds; two calls
/ with same baseSeed must produce identical paths element-by-element.

pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;7;1f%252f;0.02;0f;42);

ensA:.strategy.path.ensemble[pathCfg;5;100];
ensB:.strategy.path.ensemble[pathCfg;5;100];
.testutil.assertTrue[5=count ensA;"ensemble has numPaths paths"];
.testutil.assertTrue[5=count ensB;"second call returns same count"];

deterministic:all {(x`spot)~y`spot}'[ensA;ensB];
.testutil.assertTrue[deterministic;"two calls with same baseSeed produce identical spot vectors"];

ensC:.strategy.path.ensemble[pathCfg;5;200];
different:not (first[ensA]`spot)~first[ensC]`spot;
.testutil.assertTrue[different;"different baseSeed yields different paths"];

p0:ensA 0;
p1:ensA 1;
.testutil.assertTrue[not (p0`spot)~p1`spot;"paths within ensemble differ (seeds 100 vs 101)"];

zeroResult:.[.strategy.path.ensemble;(pathCfg;0;100);{`ERROR}];
.testutil.assertTrue[zeroResult~`ERROR;"numPaths <= 0 rejected"];

-1 "PASS test_strategy_path_ensemble: count=",string[count ensA],", deterministic=",string deterministic;
