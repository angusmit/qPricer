\l lib/init.q
/ chiReversion: known-answer hysteresis (enter long when z<-entryZ, short when
/ z>entryZ, exit inside exitZ, else hold) + independent-revaluation accounting.
n:8;
dates:2020.01.01+til n;
chiZ:-2.0 -1.5 -0.1 0.2 1.6 2.0 0.5 -1.2;
ret:0.01 -0.005 0.02 0.01 -0.01 0.015 -0.02 0.005;
p:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`volTargetScale`convenienceYield`curveSlopeCarry`chiZ`isTrain!(
    til n; dates; 60f+sums ret; n#0Nf; n#0.02; n#0f; 60f+sums ret; n#`OK; ret; 60f+sums ret; n#1f; n#0f; n#0f; chiZ; n#0b);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;60f;0.25;1f);
m:.model.createBlackScholesModel[]; fc:()!();

/ entryZ=1.0, exitZ=0.3, zero cost. Expected: -2(long),-1.5(long),-0.1(exit 0),
/ 0.2(0),1.6(short),2.0(short),0.5(hold short),-1.2(long).
cfg0:@[.strategy.defaultConfig `chiReversion;`txnCostRate;:;0f];
res0:.strategy.run[`chiReversion;trade;p;m;fc;cfg0];
.testutil.assertTrue[(1 1 0 0 -1 -1 -1 1f)~res0`rawTarget;"hysteresis enter/exit/hold logic"];

/ Accounting: independent vectorised P&L == engine stepPnl (cost on).
txnRate:0.0005;
cfg1:@[.strategy.defaultConfig `chiReversion;`txnCostRate;:;txnRate];
res1:.strategy.run[`chiReversion;trade;p;m;fc;cfg1];
positions:res1`position;
prevPos:0f,-1_positions;
turnover:(abs positions-prevPos)*txnRate;
indepStepPnl:(prevPos*ret)-turnover;
.testutil.assertTrue[1e-12>max abs indepStepPnl-res1`stepPnl;"vectorised P&L == engine stepPnl"];
.testutil.assertTrue[1e-12>max abs (res1`cumulativePnl)-sums res1`stepPnl;"cumulativePnl == sums stepPnl"];

-1 "PASS test_strategy_chi_reversion";
