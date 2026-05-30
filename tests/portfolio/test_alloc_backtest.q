\l core/init.q
/ Walk-forward OOS backtest + compare + causality for .alloc (synthetic, deterministic).
/ Three strategy return series of length 40 (different cycle lengths -> non-degenerate cov).
a:40#0.01 -0.02 0.015 -0.01 0.02 0.005;
b:40#-0.01 0.02 0.01 -0.015;
c:40#0.005 -0.005 0.02 -0.02 0.01;
r:(a;b;c);
splitCfg:`scheme`trainSpan`testSpan`maxSplits!(`rolling;15;7;4);
methods:`equalWeight`inverseVol`minVariance`riskParity`maxSharpe;
cmp:.alloc.compare[r;methods;.alloc.defaultConfig[];splitCfg];

.testutil.assertTrue[(count methods)=count cmp;"compare: one row per method"];
.testutil.assertTrue[(`method`oosSharpe`oosAnnReturn`oosAnnVol`oosMaxDrawdown`oosHitRate`avgTurnover`nSplits`oosLen`avgWeights)~cols cmp;"compare schema"];
.testutil.assertTrue[(cmp`oosSharpe)~desc cmp`oosSharpe;"ranked descending by oosSharpe"];
.testutil.assertTrue[`equalWeight in cmp`method;"equalWeight baseline present"];
.testutil.assertTrue[all (cmp`oosLen)=first cmp`oosLen;"same OOS length across methods"];
.testutil.assertTrue[all 1<cmp`nSplits;"more than one split"];

/ equalWeight OOS series == the simple mean of the strategies over the OOS windows.
ew:.alloc.backtest[r;`equalWeight;.alloc.defaultConfig[];splitCfg];
.testutil.assertTrue[not null ew`oosSharpe;"equalWeight oosSharpe finite"];
/ realized avg weights for equalWeight are exactly 1/N.
.testutil.assertTrue[1e-12>max abs (ew`avgWeights)-(count r)#1f%count r;"equalWeight realized avg weights == 1/N"];

/ CAUSALITY: weights computed on a train slice do not change when FUTURE columns are
/ appended (index-based splits, train-only estimation -> no look-ahead).
trainCols:til 15;
rExt:r,'(3 5)#0.003 -0.004 0.006 -0.002 0.001 0.004 -0.003 0.002 0.005 -0.001 0.003 -0.002 0.004 -0.005 0.001;
wShort:.alloc.weights[r[;trainCols];`minVariance;.alloc.defaultConfig[]];
wExt:.alloc.weights[rExt[;trainCols];`minVariance;.alloc.defaultConfig[]];
.testutil.assertTrue[wShort~wExt;"causality: train-slice weights unchanged when future data appended"];
/ and the split definitions are stable when the window grows (reused commodityBT.__splits).
s40:.strategy.commodityBT.__splits[`rolling;40;15;7;4];
s60:.strategy.commodityBT.__splits[`rolling;60;15;7;4];
.testutil.assertTrue[s40~(count s40)#s60;"causality: existing splits unchanged when data appended"];

-1 "PASS test_alloc_backtest: OOS compare ranks, equalWeight baseline, causal weights";
