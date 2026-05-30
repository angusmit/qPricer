\l core/init.q
/ Legacy-equivalence round-trip: a synthetic commodity-strategy run under the
/ DEFAULT execution config must reproduce the legacy flat-cost booking exactly -
/ stepPnl = prevPos*ret - txnRate*|target-prevPos|, full fill (filled==target),
/ zero slippage / fixed. This proves the .exec.fill wiring is faithful (the
/ byte-identical crux of migration step 4).
n:12;
dates:2020.01.01+til n;
ret:0.01 -0.005 0.02 0.01 -0.01 0.015 -0.02 0.005 0.008 -0.012 0.004 0.01;
chiZ:-2.0 -1.5 -0.1 0.2 1.6 2.0 0.5 -1.2 -0.3 0.4 1.1 -1.5;
isT:n#0b; isT[til 6]:1b;
p:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`volTargetScale`chiZ`isTrain!(
    til n; dates; 60f+sums ret; n#0Nf; n#0.02; n#0f; 60f+sums ret; n#`OK; ret; 60f+sums ret; n#1f; chiZ; isT);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;60f;0.25;1f);
m:.model.createBlackScholesModel[]; fc:()!();

cfg:.strategy.defaultConfig `chiReversion;
txnRate:cfg`txnCostRate;
res:.strategy.run[`chiReversion;trade;p;m;fc;cfg];

/ Independent legacy recompute from the result columns.
positions:res`position;             / held (filled) position after each step
targets:res`targetPosition;         / desired target each step
rets:res`frontReturn;
prevPos:0f,-1_positions;            / position held going into each step
expPositionPnl:prevPos*rets;
expPositionPnl[0]:0f;               / entry step books no return PnL
expCost:txnRate*abs targets-prevPos;
expStepPnl:expPositionPnl-expCost;

.testutil.assertTrue[1e-12>max abs expStepPnl-res`stepPnl;"stepPnl == independent legacy recompute"];
.testutil.assertTrue[1e-12>max abs (res`turnoverCost)-expCost;"turnoverCost == txnRate*|target-prevPos|"];
.testutil.assertTrue[(res`position)~res`targetPosition;"default config fills fully (filled == target)"];
.testutil.assertTrue[1e-12>max abs res`slippageCost;"default config has zero slippage"];
.testutil.assertTrue[1e-12>max abs res`fixedCost;"default config has zero fixed cost"];
.testutil.assertTrue[1e-12>max abs (res`proportionalCost)-res`turnoverCost;"totalCost == proportionalCost under default"];
.testutil.assertTrue[(res`filledTurnover)~res`targetTurnover;"filled turnover == target turnover under default"];

-1 "PASS test_exec_commodity_bt_legacy: default exec config reproduces legacy cost/PnL byte-identically";
