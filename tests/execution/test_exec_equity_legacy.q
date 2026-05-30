\l core/init.q
/ Equity legacy-equivalence (Step 4b): the hedge leg now routes through .exec.fill,
/ and the DEFAULT (frictionless) config reproduces the legacy flat hedge cost exactly.
/ (a) default exec == an explicit frictionless exec sub-config (byte-identical run);
/ (b) with a nonzero txnCostRate the per-step hedge cost == legacy |hedgeTrade|*spot*rate.
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0.1;0.25;9;1f%252f;0.05;0f;7);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `EQL;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
rate:0.0005;
base:@[.strategy.defaultConfig `gammaScalp;(`rebalanceMode;`rebalanceInterval;`stepYears;`txnCostRate);:;(`interval;1;1f%252f;rate)];

rDefault:(.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;base])`result;
frictionless:base,(enlist `exec)!enlist `slippageBps`fixedPerTrade`participationCap!(0f;0f;0n);
rFr:(.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;frictionless])`result;

.testutil.assertTrue[(rDefault`stepPnl)~rFr`stepPnl;"default exec == explicit frictionless exec (byte-identical stepPnl)"];
.testutil.assertTrue[(rDefault`txnCost)~rFr`txnCost;"txnCost identical under explicit frictionless"];
.testutil.assertTrue[(rDefault`hedgePosition)~rFr`hedgePosition;"hedge position identical"];

/ legacy proportional cost: per-step hedge txnCost == |hedgeTrade| * spot * rate.
expCost:(abs (rDefault`hedgeTrade)*rDefault`spot)*rate;
.testutil.assertTrue[1e-12>max abs (rDefault`txnCost)-expCost;"hedge txnCost == legacy |hedgeTrade|*spot*rate"];

-1 "PASS test_exec_equity_legacy: wired hedge reproduces legacy flat cost byte-identically";
