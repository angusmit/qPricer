\l core/init.q
/ Regime-conditional breakdown: per-bucket performance == hand-computed split, and the
/ (blended) row == the all-data __perf (the backtest headline). Synthetic, no HDB.
n:20;
dts:2020.01.02+til n;
pnl:0.01 -0.02 0.03 -0.01 0.02 0.01 -0.03 0.02 0.01 -0.01 0.04 -0.02 0.01 0.03 -0.01 0.02 -0.04 0.01 0.02 -0.01;
/ alternate the regime bucket so each has several days.
bucket:n#`backwardation`contango;
pnlByDate:([] date:dts; pnl:pnl);
regimeLabels:([] date:dts; curveState:bucket; volState:n#`normal);

br:.regime.breakdown[pnlByDate;regimeLabels;`curveState];
.testutil.assertTrue[(`bucket`nDays`totalPnl`annualReturn`annualVol`sharpe`maxDrawdown`hitRate)~cols br;"breakdown schema"];
.testutil.assertTrue[(`$"(blended)") in br`bucket;"a (blended) row is present"];
.testutil.assertTrue[(2+1)=count br;"two buckets + blended"];

annDays:.cfg.regime`annualizationDays;
/ blended row reproduces the all-data __perf exactly.
blendRow:first select from br where bucket=`$"(blended)";
expBlend:.strategy.commodityBT.__perf[pnl;annDays;1f];
.testutil.assertNear[blendRow`sharpe;expBlend`sharpe;1e-12;"(blended) sharpe == all-data __perf sharpe"];
.testutil.assertNear[blendRow`totalPnl;expBlend`totalPnl;1e-12;"(blended) totalPnl == all-data __perf"];
.testutil.assertTrue[n=blendRow`nDays;"(blended) covers all days"];

/ each bucket row == the hand split of that bucket's pnl.
bwPnl:pnl where bucket=`backwardation;
coPnl:pnl where bucket=`contango;
bwRow:first select from br where bucket=`backwardation;
coRow:first select from br where bucket=`contango;
expBw:.strategy.commodityBT.__perf[bwPnl;annDays;1f];
expCo:.strategy.commodityBT.__perf[coPnl;annDays;1f];
.testutil.assertNear[bwRow`sharpe;expBw`sharpe;1e-12;"backwardation bucket sharpe == hand split"];
.testutil.assertNear[bwRow`totalPnl;expBw`totalPnl;1e-12;"backwardation bucket totalPnl == hand split"];
.testutil.assertNear[coRow`sharpe;expCo`sharpe;1e-12;"contango bucket sharpe == hand split"];
.testutil.assertTrue[(bwRow`nDays)=count bwPnl;"backwardation nDays == split count"];
/ the buckets partition the days: bucket nDays sum to blended nDays.
.testutil.assertTrue[(blendRow`nDays)=(bwRow`nDays)+coRow`nDays;"bucket day counts partition the blended total"];

-1 "PASS test_regime_breakdown: per-bucket __perf == hand split; (blended) == backtest headline";
