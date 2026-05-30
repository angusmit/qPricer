\l core/init.q
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.30;5;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `LV_G;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
baseCfg:@[.strategy.defaultConfig `longVol;`stepYears;:;1f%252f];
closedCfg:@[baseCfg;(`forecastVol;`entryMargin);:;(0.10;0.02)];
closedBundle:.strategy.runAndSummarize[`longVol;trade;pathTbl;bsModel;fdmCfg;closedCfg];
.testutil.assertTrue[not closedBundle[`summary;`gateOpen];"forecast < implied -> closed"];
.testutil.assertTrue[`flat=closedBundle[`summary;`status];"flat status"];
.testutil.assertTrue[0f=closedBundle[`summary;`totalPnl];"totalPnl=0 when flat"];
openCfg:@[baseCfg;(`forecastVol;`entryMargin);:;(0.50;0.02)];
openBundle:.strategy.runAndSummarize[`longVol;trade;pathTbl;bsModel;fdmCfg;openCfg];
.testutil.assertTrue[openBundle[`summary;`gateOpen];"forecast 0.50 > implied 0.30 + 0.02 -> open"];
.testutil.assertTrue[(openBundle[`summary;`premiumPaid])>0f;"premium > 0 when open"];
-1 "PASS test_strategy_long_vol_gate: closedTotalPnl=",string[closedBundle[`summary;`totalPnl]];
