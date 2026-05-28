\l lib/init.q
/ performanceByStrategy aggregates per-strategy mean/std/percentiles. Cross-check
/ against an independent recomputation on a synthetic ensembleSummary.

synthRows:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`A;0;1.0;0.5;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`A;1;2.0;0.3;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`A;2;-1.0;0.7;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`A;3;3.0;0.1;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`A;4;0.5;0.6;1;`OK;"");
synthSummary:(,/) enlist each synthRows;

perfTbl:.strategy.portfolio.performanceByStrategy synthSummary;

requiredCols:`strategyName`pathCount`meanPnl`stdPnl`minPnl`maxPnl`p05Pnl`p95Pnl`sharpeLike`meanMaxDrawdown`winRate;
.testutil.assertTableColumns[perfTbl;requiredCols;"performance schema"];
.testutil.assertTrue[1=count perfTbl;"one strategy in summary"];

rowA:perfTbl 0;
.testutil.assertTrue[5=rowA`pathCount;"pathCount = 5"];

pnls:1 2 -1 3 0.5;
expectedMean:avg pnls;
expectedStd:dev pnls;
expectedMin:min pnls;
expectedMax:max pnls;
expectedWinRate:`float$(sum pnls>0f)%count pnls;
expectedMeanDD:avg 0.5 0.3 0.7 0.1 0.6;

.testutil.assertNear[rowA`meanPnl;expectedMean;1e-12;"meanPnl matches independent avg"];
.testutil.assertNear[rowA`stdPnl;expectedStd;1e-12;"stdPnl matches independent dev"];
.testutil.assertNear[rowA`minPnl;expectedMin;1e-12;"minPnl matches"];
.testutil.assertNear[rowA`maxPnl;expectedMax;1e-12;"maxPnl matches"];
.testutil.assertNear[rowA`winRate;expectedWinRate;1e-12;"winRate matches"];
.testutil.assertNear[rowA`meanMaxDrawdown;expectedMeanDD;1e-12;"meanMaxDrawdown matches"];

sortedPnls:asc pnls;
expectedP05:sortedPnls 0;
expectedP95:sortedPnls 3;
.testutil.assertNear[rowA`p05Pnl;expectedP05;1e-12;"p05Pnl independent"];
.testutil.assertNear[rowA`p95Pnl;expectedP95;1e-12;"p95Pnl independent"];

expectedSharpe:expectedMean%expectedStd;
.testutil.assertNear[rowA`sharpeLike;expectedSharpe;1e-12;"sharpeLike = mean / std"];

-1 "PASS test_strategy_portfolio_performance: meanPnl=",string[rowA`meanPnl],", expected=",string[expectedMean];
