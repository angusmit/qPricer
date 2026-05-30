\l lib/init.q
/ B2 curveRelativeValue: known-answer (long the cheap-minus-rich convergence
/ spread when the residual gap exceeds minGap) + accounting on the rv-spread
/ return column. Also verifies the rvSpreadReturn precompute on a small curve.
n:8; dates:2020.01.01+til n;
ret:0.01 -0.005 0.02 0.01 -0.01 0.015 -0.02 0.005;
rvRet:0.001 0.0008 -0.0005 0.0012 0.0003 -0.0007 0.0009 0.0004;
rvSig:0.02 0.03 0.01 0.04 0.0 0.05 0.01 0.03;
p:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`volTargetScale`rvSignal`rvSpreadReturn`volTargetScaleRV`isTrain!(
    til n; dates; 60f+sums ret; n#0Nf; n#0.02; n#0f; 60f+sums ret; n#`OK; ret; 60f+sums ret; n#1f; rvSig; rvRet; n#3f; n#0b);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;60f;0.25;1f);
m:.model.createBlackScholesModel[]; fc:()!();

/ minGap default 0 -> long whenever gap>0 (flat when gap=0).
cfg0:@[.strategy.defaultConfig `curveRelativeValue;`txnCostRate;:;0f];
r0:.strategy.run[`curveRelativeValue;trade;p;m;fc;cfg0];
expectedRaw:?[rvSig>0f;1f;0f];
.testutil.assertTrue[expectedRaw~r0`rawTarget;"long the convergence spread when residual gap > minGap"];

txnRate:0.0005;
cfg1:@[.strategy.defaultConfig `curveRelativeValue;`txnCostRate;:;txnRate];
r1:.strategy.run[`curveRelativeValue;trade;p;m;fc;cfg1];
positions:r1`position; prevPos:0f,-1_positions;
indep:(prevPos*rvRet)-(abs positions-prevPos)*txnRate;
.testutil.assertTrue[1e-12>max abs indep-r1`stepPnl;"vectorised rv-spread P&L == engine stepPnl"];
.testutil.assertTrue[1e-12>max abs (r1`cumulativePnl)-sums r1`stepPnl;"cumulativePnl == sums stepPnl"];

/ rvSpreadReturn precompute: on a 3-tenor curve, cheap/rich legs are identified
/ and the spread return is finite (no crash). Use a curve with curvature.
contracts:202003 202004 202005 202006;
expiries:2020.03.20 2020.04.20 2020.05.20 2020.06.20;
mkRows:{[d;contracts;expiries]
    alive:where expiries>=d;
    tens:(`float$(expiries alive)-d)%365f;
    base:55f-2.5f*tens+0.4f*sin 25f*tens;
    ([] asofDate:(count alive)#d; tenor:tens; price:base; contractYM:contracts alive; expiry:expiries alive)};
hist:raze mkRows[;contracts;expiries] each 2020.01.20+til 12;
sig:.strategy.path.commoditySignals[hist;`trainEndDate`kalmanEstCfg!(2020.01.25;`gridSteps`refineRounds`nSweeps!(3;1;1))];
.testutil.assertTrue[all not null (sig`path)`rvSpreadReturn;"rvSpreadReturn finite on a real-shaped curve"];

-1 "PASS test_strategy_curve_relative_value";
