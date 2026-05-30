\l core/init.q
/ B3 twoTimescale: known-answer (weighted thresholded chi-reversion + xi-momentum
/ trend) + independent-revaluation accounting.
n:8; dates:2020.01.01+til n;
ret:0.01 -0.005 0.02 0.01 -0.01 0.015 -0.02 0.005;
chiZ:-2.0 -1.5 -0.1 0.2 1.6 2.0 0.5 -1.2;
xiMom:0.5 0.3 -0.2 -0.4 0.1 0.2 -0.3 0.4;
p:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`volTargetScale`chiZ`xiMomentum`isTrain!(
    til n; dates; 60f+sums ret; n#0Nf; n#0.02; n#0f; 60f+sums ret; n#`OK; ret; 60f+sums ret; n#1f; chiZ; xiMom; n#0b);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;60f;0.25;1f);
m:.model.createBlackScholesModel[]; fc:()!();

/ revertWeight=trendWeight=0.5, entryZ=1.0 -> combined = 0.5*revert + 0.5*trend.
cfg0:@[.strategy.defaultConfig `twoTimescale;`txnCostRate;:;0f];
r0:.strategy.run[`twoTimescale;trade;p;m;fc;cfg0];
revertRaw:?[chiZ<-1.0;1f;?[chiZ>1.0;-1f;0f]];
expectedRaw:(0.5*revertRaw)+0.5*signum xiMom;
.testutil.assertTrue[expectedRaw~r0`rawTarget;"weighted revert+trend known-answer"];

txnRate:0.0005;
cfg1:@[.strategy.defaultConfig `twoTimescale;`txnCostRate;:;txnRate];
r1:.strategy.run[`twoTimescale;trade;p;m;fc;cfg1];
positions:r1`position; prevPos:0f,-1_positions;
indep:(prevPos*ret)-(abs positions-prevPos)*txnRate;
.testutil.assertTrue[1e-12>max abs indep-r1`stepPnl;"vectorised P&L == engine stepPnl"];
.testutil.assertTrue[1e-12>max abs (r1`cumulativePnl)-sums r1`stepPnl;"cumulativePnl == sums stepPnl"];

-1 "PASS test_strategy_two_timescale";
