\l lib/init.q
/ convenienceYieldCarry: known-answer position sign (long when convenience yield
/ exceeds the rate = backwardation; short when below = contango) + independent-
/ revaluation accounting (vectorised position*return P&L == engine stepPnl).
n:8;
dates:2020.01.01+til n;
cy:0.10 0.08 0.06 0.04 -0.02 -0.04 -0.06 -0.08;
ret:0.01 -0.005 0.02 0.01 -0.01 0.015 -0.02 0.005;
buildPath:{[n;dates;cy;ret]
    flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`volTargetScale`convenienceYield`curveSlopeCarry`chiZ`isTrain!(
        til n; dates; 60f+sums ret; n#0Nf; n#0.02; n#0f; 60f+sums ret; n#`OK; ret; 60f+sums ret; n#1f; cy; cy; n#0f; n#0b)};
p:buildPath[n;dates;cy;ret];
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;60f;0.25;1f);
m:.model.createBlackScholesModel[]; fc:()!();

/ Known-answer: zero cost -> long while cy>rate, short while cy<rate.
cfg0:@[.strategy.defaultConfig `convenienceYieldCarry;`txnCostRate;:;0f];
res0:.strategy.run[`convenienceYieldCarry;trade;p;m;fc;cfg0];
.testutil.assertTrue[(1 1 1 1 -1 -1 -1 -1f)~res0`rawTarget;"long backwardation / short contango"];

/ Accounting: independent vectorised recompute == engine stepPnl (with cost on).
txnRate:0.0005;
cfg1:@[.strategy.defaultConfig `convenienceYieldCarry;`txnCostRate;:;txnRate];
res1:.strategy.run[`convenienceYieldCarry;trade;p;m;fc;cfg1];
positions:res1`position;
prevPos:0f,-1_positions;
turnover:(abs positions-prevPos)*txnRate;
indepStepPnl:(prevPos*ret)-turnover;
.testutil.assertTrue[1e-12>max abs indepStepPnl-res1`stepPnl;"vectorised P&L == engine stepPnl"];
/ Portfolio-value identity within the result: cumulativePnl == running sum of stepPnl.
.testutil.assertTrue[1e-12>max abs (res1`cumulativePnl)-sums res1`stepPnl;"cumulativePnl == sums stepPnl (PV=cash identity)"];

-1 "PASS test_strategy_convenience_carry";
