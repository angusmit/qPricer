\l core/init.q
/ strategyCorrelation: symmetric, unit diagonal, perfect self-correlation = 1.

synthRows:();
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`X;0;1.0;0f;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`X;1;2.0;0f;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`X;2;3.0;0f;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`Y;0;2.0;0f;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`Y;1;4.0;0f;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`Y;2;6.0;0f;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`Z;0;-1.0;0f;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`Z;1;-2.0;0f;1;`OK;"");
synthRows,:enlist `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(`Z;2;-3.0;0f;1;`OK;"");
synthSummary:(,/) enlist each synthRows;

corDict:.strategy.portfolio.strategyCorrelation synthSummary;
sNames:corDict`names;
sMatrix:corDict`matrix;

.testutil.assertTrue[(asc `X`Y`Z)~asc sNames;"three strategies in correlation"];
.testutil.assertTrue[3=count sMatrix;"3x3 matrix rows"];

xIdx:sNames?`X;
yIdx:sNames?`Y;
zIdx:sNames?`Z;

.testutil.assertNear[sMatrix[xIdx;xIdx];1f;1e-12;"self-correlation X-X = 1"];
.testutil.assertNear[sMatrix[yIdx;yIdx];1f;1e-12;"self-correlation Y-Y = 1"];
.testutil.assertNear[sMatrix[zIdx;zIdx];1f;1e-12;"self-correlation Z-Z = 1"];

.testutil.assertNear[sMatrix[xIdx;yIdx];1f;1e-12;"X and Y are scaled copies (Y=2X) -> corr = 1"];
.testutil.assertNear[sMatrix[xIdx;zIdx];-1f;1e-12;"X and Z are anti-correlated (Z=-X) -> corr = -1"];

.testutil.assertNear[sMatrix[xIdx;yIdx];sMatrix[yIdx;xIdx];1e-12;"matrix symmetric (X,Y) = (Y,X)"];
.testutil.assertNear[sMatrix[xIdx;zIdx];sMatrix[zIdx;xIdx];1e-12;"matrix symmetric (X,Z) = (Z,X)"];

-1 "PASS test_strategy_portfolio_correlation: XY=",string[sMatrix[xIdx;yIdx]],", XZ=",string[sMatrix[xIdx;zIdx]];
