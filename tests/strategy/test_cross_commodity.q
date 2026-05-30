\l lib/init.q
/ Cross-commodity aggregation arithmetic vs an independent recompute on a hand
/ set of per-(commodity,split) numbers, and a stable descending rank.
clDetail:([] strategyName:`a`a`b`b; splitId:0 1 0 1;
    testSharpe:1.0 0.5 -0.5 2.0; testAnnualReturn:0.1 0.05 -0.05 0.2; testMaxDrawdown:0.02 0.03 0.04 0.01);
ngDetail:([] strategyName:`a`a`b`b; splitId:0 1 0 1;
    testSharpe:0.0 -1.0 1.5 0n; testAnnualReturn:0.0 -0.1 0.15 0.0; testMaxDrawdown:0.05 0.06 0.02 0.07);
agg:.strategy.commodityBT.crossCommodity[`CL`NG!(clDetail;ngDetail)];

/ a: sharpes across both commodities = 1.0 0.5 0.0 -1.0 ; b: -0.5 2.0 1.5 (null skipped).
aRow:first select from agg where strategyName=`a;
bRow:first select from agg where strategyName=`b;
.testutil.assertNear[aRow`meanSharpe;avg 1.0 0.5 0.0 -1.0;1e-12;"a cross-commodity mean Sharpe"];
.testutil.assertNear[aRow`sharpeStd;dev 1.0 0.5 0.0 -1.0;1e-12;"a dispersion"];
.testutil.assertTrue[4=aRow`nCells;"a has 4 (commodity x split) cells"];
.testutil.assertNear[aRow`fractionPositive;2%4;1e-12;"a fraction of positive cells (2 of 4)"];
.testutil.assertNear[bRow`meanSharpe;avg -0.5 2.0 1.5;1e-12;"b mean Sharpe skips the null cell"];
.testutil.assertTrue[3=bRow`nTraded;"b traded in 3 of 4 cells"];
.testutil.assertNear[bRow`fractionPositive;2%4;1e-12;"b fraction positive (null is not positive)"];

/ Ranked descending by meanSharpe (b mean 1.0 > a mean 0.125).
.testutil.assertTrue[(agg`meanSharpe)~desc agg`meanSharpe;"ranked descending by mean Sharpe"];
.testutil.assertTrue[`b~first agg`strategyName;"b ranks above a"];

-1 "PASS test_cross_commodity";
