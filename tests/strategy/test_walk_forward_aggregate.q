\l lib/init.q
/ Walk-forward aggregation arithmetic vs an independent recompute on a hand set
/ of per-split numbers (avg / dev skip nulls = strategies that never traded).
detail:([] strategyName:`a`a`a`b`b`b;
    splitId:0 1 2 0 1 2;
    testSharpe:1.0 -0.5 0.6 2.0 0n 1.0;
    testAnnualReturn:0.10 -0.05 0.06 0.20 0.0 0.10;
    testMaxDrawdown:0.02 0.03 0.01 0.01 0.04 0.02);
agg:.strategy.commodityBT.__aggregateSplits detail;

aRow:first select from agg where strategyName=`a;
.testutil.assertNear[aRow`meanSharpe;avg 1.0 -0.5 0.6;1e-12;"a meanSharpe == independent avg"];
.testutil.assertNear[aRow`sharpeStd;dev 1.0 -0.5 0.6;1e-12;"a dispersion == independent dev"];
.testutil.assertTrue[2=aRow`splitsPositive;"a splits-positive count (2 of 3)"];
.testutil.assertTrue[3=aRow`nTraded;"a traded in all 3 splits"];

bRow:first select from agg where strategyName=`b;
.testutil.assertNear[bRow`meanSharpe;avg 2.0 1.0;1e-12;"b meanSharpe skips the null split"];
.testutil.assertTrue[2=bRow`nTraded;"b traded in 2 of 3 splits (one null)"];
.testutil.assertTrue[2=bRow`splitsPositive;"b splits-positive count (null is not positive)"];
.testutil.assertTrue[3=bRow`nSplits;"b appears in 3 splits total"];

-1 "PASS test_walk_forward_aggregate";
