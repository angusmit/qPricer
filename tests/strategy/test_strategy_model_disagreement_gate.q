\l core/init.q
/ Gate logic: open when |disagreement| > threshold; closed otherwise. Direction = auto
/ follows sign; forced long/short overrides.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;5;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `MD_G;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);

baseCfg:.strategy.defaultConfig `modelDisagreement;
baseCfg:@[baseCfg;`stepYears;:;1f%252f];

closedCfg:@[baseCfg;(`modelBVolBump;`disagreementThreshold);:;(0.01;0.5)];
closedBundle:.strategy.runAndSummarize[`modelDisagreement;trade;pathTbl;bsModel;fdmCfg;closedCfg];
.testutil.assertTrue[not closedBundle[`summary;`gateOpen];"small disagreement < threshold -> closed"];
.testutil.assertTrue[`flat=closedBundle[`summary;`status];"summary flat status"];
.testutil.assertTrue[0f=closedBundle[`summary;`totalPnl];"flat totalPnl 0"];

openCfg:@[baseCfg;(`modelBVolBump;`disagreementThreshold);:;(0.10;0.05)];
openBundle:.strategy.runAndSummarize[`modelDisagreement;trade;pathTbl;bsModel;fdmCfg;openCfg];
.testutil.assertTrue[openBundle[`summary;`gateOpen];"big disagreement > threshold -> open"];
.testutil.assertTrue[1=openBundle[`summary;`tradeSide];"auto direction goes long when modelB > modelA"];

negCfg:@[openCfg;`modelBVolBump;:;-0.10];
negBundle:.strategy.runAndSummarize[`modelDisagreement;trade;pathTbl;bsModel;fdmCfg;negCfg];
.testutil.assertTrue[-1=negBundle[`summary;`tradeSide];"auto direction goes short when modelB < modelA"];

forcedCfg:@[openCfg;`direction;:;`short];
forcedBundle:.strategy.runAndSummarize[`modelDisagreement;trade;pathTbl;bsModel;fdmCfg;forcedCfg];
.testutil.assertTrue[-1=forcedBundle[`summary;`tradeSide];"forced short overrides auto"];

-1 "PASS test_strategy_model_disagreement_gate: closedTotalPnl=",string[closedBundle[`summary;`totalPnl]],", openTradeSide=",string[openBundle[`summary;`tradeSide]];
