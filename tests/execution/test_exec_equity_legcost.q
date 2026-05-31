\l core/init.q
/ Step 4c: equity OPTION-LEG costs route through .exec.fill. Legacy-equivalence per
/ form (default == legacy, byte-identical) + realism (premium-scaled option-leg
/ slippage; net<gross; deltaPV==stepPnl identity preserved). iron_condor is held
/ outright (no hedge by default) so its entire cost is the 4-leg option entry.
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;6;1f%252f;0.02;0f;29);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `LC;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
rate:0.001;
base:@[.strategy.defaultConfig `ironCondor;(`stepYears;`financingRate;`txnCostRate);:;(1f%252f;0f;rate)];

/ (a) legacy-equivalence: default exec == an explicit frictionless exec sub-config.
rDef:(.strategy.runAndSummarize[`ironCondor;trade;pathTbl;bsModel;fdmCfg;base])`result;
fr:base,(enlist `exec)!enlist `slippageBps`fixedPerTrade`participationCap!(0f;0f;0n);
rFr:(.strategy.runAndSummarize[`ironCondor;trade;pathTbl;bsModel;fdmCfg;fr])`result;
.testutil.assertTrue[(rDef`stepPnl)~rFr`stepPnl;"iron_condor leg cost: default == explicit frictionless (byte-identical)"];
.testutil.assertTrue[(rDef`txnCost)~rFr`txnCost;"iron_condor: leg txnCost identical under frictionless"];

/ (b) realism: option-leg slippage. net < gross; gross-net == slippage paid; linear in bps.
net50:base,(enlist `exec)!enlist (enlist `slippageBps)!enlist 50f;
net100:base,(enlist `exec)!enlist (enlist `slippageBps)!enlist 100f;
rN50:(.strategy.runAndSummarize[`ironCondor;trade;pathTbl;bsModel;fdmCfg;net50])`result;
rN100:(.strategy.runAndSummarize[`ironCondor;trade;pathTbl;bsModel;fdmCfg;net100])`result;
gTot:sum rDef`stepPnl; nTot:sum rN50`stepPnl;
slip50:(sum rN50`txnCost)-sum rDef`txnCost;
slip100:(sum rN100`txnCost)-sum rDef`txnCost;
.testutil.assertTrue[slip50>0f;"option-leg slippage is paid"];
.testutil.assertTrue[nTot<gTot;"net PnL < gross PnL under option-leg slippage"];
.testutil.assertTrue[1e-9>abs (gTot-nTot)-slip50;"gross - net == option-leg slippage paid"];
.testutil.assertTrue[1e-9>abs slip100-2f*slip50;"leg slippage scales linearly with bps (|premium*qty|*bps)"];

/ (c) identity under option-leg slippage: deltaPV == stepPnl (fold to states).
pathRows:0!pathTbl; firstStep:first pathRows; remaining:1_pathRows;
initS:.strategy.ironCondor.init[trade;firstStep;bsModel;fdmCfg;net50];
stepFn:.strategy.ironCondor.step[;;trade;bsModel;fdmCfg;net50];
states:enlist[initS],stepFn\[initS;remaining];
spots:pathRows`spot;
pv:{[s;sp] .strategy.__portfolioValue[s`cash;s`prevPositionValue;s`hedgePosition;sp]}'[states;spots];
dPV:1_(pv-prev pv);
.testutil.assertTrue[1e-8>max abs dPV-1_rN50`stepPnl;"deltaPV==stepPnl identity holds under option-leg slippage"];

/ (d) a DIFFERENT leg-cost form (calendar roll close/open) is also byte-identical by default.
crBase:@[.strategy.defaultConfig `calendarRoll;(`stepYears;`txnCostRate);:;(1f%252f;0.001)];
crFr:crBase,(enlist `exec)!enlist `slippageBps`fixedPerTrade`participationCap!(0f;0f;0n);
crPath:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0.05;0.2;30;1f%252f;0.02;0f;17);
rcr0:(.strategy.runAndSummarize[`calendarRoll;trade;crPath;bsModel;fdmCfg;crBase])`result;
rcr1:(.strategy.runAndSummarize[`calendarRoll;trade;crPath;bsModel;fdmCfg;crFr])`result;
.testutil.assertTrue[(rcr0`stepPnl)~rcr1`stepPnl;"calendarRoll roll-cost form: default == frictionless byte-identical"];

-1 "PASS test_exec_equity_legcost: option-leg costs route through .exec.fill; byte-identical default + premium-scaled slippage + identity preserved";
