\l core/init.q
/ Equity realism (Step 4b): gammaScalp with hedge slippage. The deltaPV==stepPnl
/ identity STILL holds (cost threaded through cash), gross-net == the slippage paid,
/ net PnL < gross, and a participation cap makes the hedge lag its target.
pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0.15;0.30;9;1f%252f;0.05;0f;13);
pathTbl:.strategy.path.fromSynthetic pathCfg;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `EQR;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
rate:0.0005; slipBps:30f;
base:@[.strategy.defaultConfig `gammaScalp;(`rebalanceMode;`rebalanceInterval;`stepYears;`txnCostRate);:;(`interval;1;1f%252f;rate)];
net:base,(enlist `exec)!enlist (enlist `slippageBps)!enlist slipBps;

/ --- identity under slippage: fold the wired step to get states, recompute deltaPV ---
pathRows:0!pathTbl; firstStep:first pathRows; remaining:1_pathRows;
initS:.strategy.gammaScalp.init[trade;firstStep;bsModel;fdmCfg;net];
stepFn:.strategy.gammaScalp.step[;;trade;bsModel;fdmCfg;net];
states:enlist[initS],stepFn\[initS;remaining];
spots:pathRows`spot;
pvFn:{[s;sp] .strategy.__portfolioValue[s`cash;s`prevPositionValue;s`hedgePosition;sp]};
pvSeries:pvFn'[states;spots];
deltaPV:1_(pvSeries-prev pvSeries);
resNet:(.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;net])`result;
stepPnls:1_resNet`stepPnl;
idResidual:max abs deltaPV-stepPnls;
.testutil.assertTrue[idResidual<1e-8;"deltaPV==stepPnl identity holds under hedge slippage"];

/ --- gross (frictionless) vs net (slippage) ---
resGross:(.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;base])`result;
grossTot:sum resGross`stepPnl; netTot:sum resNet`stepPnl;
.testutil.assertTrue[netTot<grossTot;"net PnL < gross PnL (hedge slippage is a drag)"];
/ slippage paid = net.txnCost - gross.txnCost (same trades) == |hedge notional|*bps/1e4.
bpsFrac:slipBps%1e4;
expSlip:(abs (resNet`hedgeTrade)*resNet`spot)*bpsFrac;
slipPaid:(resNet`txnCost)-resGross`txnCost;
.testutil.assertTrue[1e-10>max abs slipPaid-expSlip;"per-step slippage == |hedge notional|*bps/1e4"];
.testutil.assertTrue[1e-10>abs (grossTot-netTot)-sum slipPaid;"gross - net == total slippage paid"];

/ --- participation cap makes the hedge lag target (direct __hedgeStep, with bar volume) ---
hs:`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(0f;0f;0f;0;100f;0f);
si:`spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand`barVolume`exec!(
    100f;0f;1f;0;1f%252f;rate;0f;`interval;1;0.05;10f;(enlist `participationCap)!enlist 0.05);
hr:.strategy.__hedgeStep[hs;si];
.testutil.assertTrue[(abs hr`hedgePosition)<abs neg si`positionDelta;"participation cap: filled hedge lags target delta (residual exposure)"];
.testutil.assertTrue[(abs hr`hedgePosition)<1e-9+0.05*10f;"filled hedge bounded by cap*barVolume (in notional/spot)"];

-1 "PASS test_exec_equity_realism: slippage drag + identity preserved + cap lags target";
