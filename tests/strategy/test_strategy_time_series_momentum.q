\l core/init.q
/ B1 timeSeriesMomentum: known-answer (position = sign of trailing momentum) +
/ independent-revaluation accounting (vectorised position*return == engine stepPnl).
n:8; dates:2020.01.01+til n;
ret:0.01 -0.005 0.02 0.01 -0.01 0.015 -0.02 0.005;
mom:0.01 0.005 -0.002 -0.01 0.003 0.008 -0.004 0.006;
p:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`volTargetScale`momentum`isTrain!(
    til n; dates; 60f+sums ret; n#0Nf; n#0.02; n#0f; 60f+sums ret; n#`OK; ret; 60f+sums ret; n#1f; mom; n#0b);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;60f;0.25;1f);
m:.model.createBlackScholesModel[]; fc:()!();

cfg0:@[.strategy.defaultConfig `timeSeriesMomentum;`txnCostRate;:;0f];
r0:.strategy.run[`timeSeriesMomentum;trade;p;m;fc;cfg0];
.testutil.assertTrue[(`float$signum mom)~r0`rawTarget;"position follows the sign of trailing momentum"];

txnRate:0.0005;
cfg1:@[.strategy.defaultConfig `timeSeriesMomentum;`txnCostRate;:;txnRate];
r1:.strategy.run[`timeSeriesMomentum;trade;p;m;fc;cfg1];
positions:r1`position; prevPos:0f,-1_positions;
indep:(prevPos*ret)-(abs positions-prevPos)*txnRate;
.testutil.assertTrue[1e-12>max abs indep-r1`stepPnl;"vectorised P&L == engine stepPnl"];
.testutil.assertTrue[1e-12>max abs (r1`cumulativePnl)-sums r1`stepPnl;"cumulativePnl == sums stepPnl"];

-1 "PASS test_strategy_time_series_momentum";
