\l lib/init.q
/ B4 storageCashCarry: known-answer (long the near-far calendar spread only when
/ contango exceeds the storage cost) + accounting on the spread-return column.
n:8; dates:2020.01.01+til n;
ret:0.01 -0.005 0.02 0.01 -0.01 0.015 -0.02 0.005;
nfRet:0.002 -0.001 0.003 0.001 -0.002 0.0015 -0.001 0.0008;
ccarry:0.05 0.03 -0.02 -0.04 -0.06 0.02 0.01 -0.03;
p:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`volTargetScale`curveSlopeCarry`nearFarSpreadReturn`volTargetScaleNearFar`isTrain!(
    til n; dates; 60f+sums ret; n#0Nf; n#0.02; n#0f; 60f+sums ret; n#`OK; ret; 60f+sums ret; n#1f; ccarry; nfRet; n#2f; n#0b);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;60f;0.25;1f);
m:.model.createBlackScholesModel[]; fc:()!();

/ storageCostRate default 0.01; contango = -curveSlopeCarry; long when contango>0.01.
cfg0:@[.strategy.defaultConfig `storageCashCarry;`txnCostRate;:;0f];
r0:.strategy.run[`storageCashCarry;trade;p;m;fc;cfg0];
expectedRaw:?[(neg ccarry)>0.01;1f;0f];
.testutil.assertTrue[expectedRaw~r0`rawTarget;"long calendar spread only when contango > storage cost"];

/ Accounting on the spread-return column (returnColumn=nearFarSpreadReturn).
txnRate:0.0005;
cfg1:@[.strategy.defaultConfig `storageCashCarry;`txnCostRate;:;txnRate];
r1:.strategy.run[`storageCashCarry;trade;p;m;fc;cfg1];
positions:r1`position; prevPos:0f,-1_positions;
indep:(prevPos*nfRet)-(abs positions-prevPos)*txnRate;
.testutil.assertTrue[1e-12>max abs indep-r1`stepPnl;"vectorised spread P&L == engine stepPnl"];
.testutil.assertTrue[1e-12>max abs (r1`cumulativePnl)-sums r1`stepPnl;"cumulativePnl == sums stepPnl"];

-1 "PASS test_strategy_storage_cash_carry";
