\l core/init.q
/ ============================================================================
/ execution_gross_vs_net_equity.q - example (NOT a test). The equity counterpart
/ of execution_gross_vs_net.q: how much of a delta-hedged equity strategy's edge
/ survives HEDGE-rebalancing slippage (v0.61, Step 4b). Runs gamma-scalping on a
/ synthetic path at rising slippage levels and reports GROSS (frictionless) vs NET
/ total PnL + a Sharpe-like ratio. Synthetic -> always runs, no real data needed.
/ ============================================================================
pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0.08;0.30;60;1f%252f;0.05;0f;42);
pathTbl:.strategy.path.fromSynthetic pathCfg;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `GS;`X;`equityOption;`european;`call;100f;0.5;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;60;100;0f;300f;`linear;1b;1b);
base:@[.strategy.defaultConfig `gammaScalp;(`rebalanceMode;`rebalanceInterval;`stepYears;`txnCostRate);:;(`interval;1;1f%252f;0.0005)];

statsAt:{[base;trade;path;m;fdm;bps]
    cfg:base,(enlist `exec)!enlist (enlist `slippageBps)!enlist bps;
    s:(.strategy.runAndSummarize[`gammaScalp;trade;path;m;fdm;cfg])`summary;
    sharpeLike:$[0f<s`stepPnlVol; (sqrt 252f)*(s`meanStepPnl)%s`stepPnlVol; 0n];
    `totalPnl`sharpeLike!((s`totalPnl);sharpeLike)}[base;trade;pathTbl;bsModel;fdmCfg];

bpsLevels:0 5 10 20 50f;
rows:statsAt each bpsLevels;
report:([] slippageBps:bpsLevels; totalPnl:rows[;`totalPnl]; sharpeLike:rows[;`sharpeLike]);
-1 "gamma-scalp GROSS vs NET by hedge slippage (synthetic path; slippageBps=0 is GROSS):";
show report;
-1 "";
-1 "gammaScalp gross totalPnl=",(string rows[0]`totalPnl),"  ->  net @20bps=",string rows[3]`totalPnl;
-1 "(Delta-hedge turnover is frequent, so hedge slippage compounds across rebalances - the honest net-of-execution view.)";
exit 0;
