\l core/init.q
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;5;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CT_M;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
collarCfg:@[.strategy.defaultConfig `collarTailHedge;`stepYears;:;1f%252f];
tailCfg:@[collarCfg;(`mode;`premiumBudgetPct);:;(`tailHedge;0.02)];
collarBundle:.strategy.runAndSummarize[`collarTailHedge;trade;pathTbl;bsModel;fdmCfg;collarCfg];
tailBundle:.strategy.runAndSummarize[`collarTailHedge;trade;pathTbl;bsModel;fdmCfg;tailCfg];
.testutil.assertTrue[`collar=collarBundle[`summary;`mode];"collar mode"];
.testutil.assertTrue[`tailHedge=tailBundle[`summary;`mode];"tailHedge mode"];
collarFirst:(collarBundle`result)0;
tailFirst:(tailBundle`result)0;
.testutil.assertTrue[(collarFirst`underlyingValue)>0f;"collar has underlying"];
.testutil.assertNear[tailFirst`underlyingValue;0f;1e-12;"tailHedge no underlying"];
.testutil.assertNear[tailFirst`callValue;0f;1e-12;"tailHedge no call"];
.testutil.assertTrue[(tailFirst`putValue)>0f;"tailHedge has puts"];
.testutil.assertNear[tailBundle[`summary;`budgetUsed];0.02;1e-3;"budgetUsed near 0.02 of spot"];
-1 "PASS test_strategy_collar_tail_hedge_mode: collarUnderlying=",string[collarFirst`underlyingValue],", tailPutValue=",string[tailFirst`putValue];
