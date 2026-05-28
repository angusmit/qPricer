\l lib/init.q
/ riskReversal entry gate: deviation = skewSlope - fairSkew. Open when |deviation| > margin.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;5;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `RR_G;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);

openCfg:.strategy.defaultConfig `riskReversal;
openCfg:@[openCfg;(`skewSlope;`fairSkew;`skewMargin;`stepYears);:;(-0.6;-0.2;0.05;1f%252f)];
openBundle:.strategy.runAndSummarize[`riskReversal;trade;pathTbl;bsModel;fdmCfg;openCfg];
.testutil.assertTrue[openBundle[`summary;`gateOpen];"deviation 0.4 > margin 0.05 -> gate open"];

closedCfg:@[openCfg;(`skewSlope;`fairSkew;`skewMargin);:;(-0.21;-0.2;0.05)];
closedBundle:.strategy.runAndSummarize[`riskReversal;trade;pathTbl;bsModel;fdmCfg;closedCfg];
closedSummary:closedBundle`summary;
closedResult:closedBundle`result;
.testutil.assertTrue[not closedSummary`gateOpen;"deviation 0.01 < margin 0.05 -> gate closed"];
.testutil.assertTrue[all (closedResult`status)=`flat;"all rows flat when gate closed"];
.testutil.assertTrue[0f=closedSummary`totalPnl;"totalPnl = 0 when flat"];
.testutil.assertTrue[`flat=closedSummary`status;"summary flat status"];

forcedCallCfg:@[openCfg;`riskReversalDirection;:;`longCallWing];
forcedCallBundle:.strategy.runAndSummarize[`riskReversal;trade;pathTbl;bsModel;fdmCfg;forcedCallCfg];
forcedFirst:(forcedCallBundle`result)0;
.testutil.assertTrue[forcedFirst[`callVol]<forcedFirst`putVol;"longCallWing forced direction valid"];

-1 "PASS test_strategy_risk_reversal_gate: openGate=",string[openBundle[`summary;`gateOpen]],", closedGate=",string[closedSummary`gateOpen];
