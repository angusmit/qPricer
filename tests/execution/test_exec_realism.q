\l core/init.q
/ Realism check: with nonzero slippage + a binding participation cap, the cost
/ components recompute independently, filled turnover <= target turnover (the cap
/ binds), slippage is paid, and NET PnL < GROSS PnL. Synthetic path carries a
/ `volume` column so the cap engages.
n:14;
dates:2020.01.01+til n;
ret:0.01 -0.02 0.03 -0.01 0.02 -0.03 0.01 0.02 -0.02 0.03 -0.01 0.02 -0.02 0.01;
chiZ:-2.0 2.0 -1.8 1.9 -2.0 1.7 -1.5 1.6 -1.9 1.8 -1.6 1.5 -1.7 2.0;
isT:n#0b; isT[til 7]:1b;
vol:n#50f;
p:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`volTargetScale`chiZ`isTrain`volume!(
    til n; dates; 60f+sums ret; n#0Nf; n#0.02; n#0f; 60f+sums ret; n#`OK; ret; 60f+sums ret; n#1f; chiZ; isT; vol);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;60f;0.25;1f);
m:.model.createBlackScholesModel[]; fc:()!();

slipBps:50f; cap:0.01;
cfg:(.strategy.defaultConfig `chiReversion),(enlist `exec)!enlist `slippageBps`participationCap!(slipBps;cap);
res:.strategy.run[`chiReversion;trade;p;m;fc;cfg];

/ cost components add up per step.
.testutil.assertTrue[1e-12>max abs (res`turnoverCost)-(res`proportionalCost)+(res`slippageCost)+res`fixedCost;"totalCost == prop + slip + fixed per step"];
/ slippage is actually paid.
.testutil.assertTrue[0f<sum res`slippageCost;"slippage cost is paid under nonzero bps"];
/ independent slippage recompute: |filled| * refPrice * bps/1e4.
expSlip:(res`filledTurnover)*(res`frontPrice)*slipBps%1e4;
.testutil.assertTrue[1e-12>max abs (res`slippageCost)-expSlip;"slippageCost == |filled|*refPrice*bps/1e4"];
/ participation cap: filled turnover <= target turnover, and the cap binds somewhere.
turnGap:(res`targetTurnover)-res`filledTurnover;
.testutil.assertTrue[all turnGap>=neg 1e-12;"filled turnover <= target turnover"];
.testutil.assertTrue[any turnGap>1e-9;"the participation cap binds (partial fills occur)"];
/ independent cap recompute: filled == min(target, cap*volume).
expFilled:(res`targetTurnover)&cap*vol;
.testutil.assertTrue[1e-12>max abs (res`filledTurnover)-expFilled;"filled == min(target, cap*volume)"];
/ NET < GROSS: execution costs strictly reduce total PnL.
netTot:sum res`stepPnl;
grossTot:sum (res`stepPnl)+res`turnoverCost;
.testutil.assertTrue[netTot<grossTot;"net PnL < gross PnL (costs reduce PnL)"];
.testutil.assertTrue[1e-12>abs (grossTot-netTot)-sum res`turnoverCost;"gross - net == total execution cost"];

-1 "PASS test_exec_realism: slip+cap -> components recompute, cap binds, net<gross";
