\l lib/init.q
/ dashboard composes performance + correlation + bookAggregate. bookAggregate is the
/ per-path SUM of per-strategy P&L (= total book P&L per path); mean/std across paths.

synthRows:();
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`A;0;1.0;0.1;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`A;1;2.0;0.2;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`A;2;3.0;0.3;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`B;0;0.5;0.05;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`B;1;1.0;0.10;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`B;2;1.5;0.15;1;`OK;"");
synthSummary:(,/) enlist each synthRows;

dash:.strategy.portfolio.dashboard synthSummary;

expectedKeys:`performanceByStrategy`strategyCorrelation`bookAggregate`dashboardStatus`dashboardMessage;
.testutil.assertTrue[all expectedKeys in key dash;"dashboard has expected keys"];
.testutil.assertTrue[`OK=dash`dashboardStatus;"status OK"];

bookAgg:dash`bookAggregate;
.testutil.assertTrue[3=bookAgg`pathCount;"bookAggregate pathCount = 3"];

bookPerPath:1.5 3.0 4.5;
.testutil.assertNear[bookAgg`meanBookPnl;avg bookPerPath;1e-12;"meanBookPnl = avg of summed per-path PnL"];
.testutil.assertNear[bookAgg`stdBookPnl;dev bookPerPath;1e-12;"stdBookPnl = dev of summed per-path PnL"];

perfTbl:dash`performanceByStrategy;
.testutil.assertTrue[2=count perfTbl;"performance has 2 strategy rows"];

corDict:dash`strategyCorrelation;
.testutil.assertTrue[(asc `A`B)~asc corDict`names;"correlation names match"];

emptyDash:.strategy.portfolio.dashboard ();
.testutil.assertTrue[`ERROR=emptyDash`dashboardStatus;"empty ensembleSummary -> ERROR status"];

-1 "PASS test_strategy_portfolio_dashboard: meanBookPnl=",string[bookAgg`meanBookPnl],", stdBookPnl=",string[bookAgg`stdBookPnl];
