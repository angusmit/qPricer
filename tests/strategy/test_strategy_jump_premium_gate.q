\l lib/init.q
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.20;5;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `JP_G;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
baseCfg:@[.strategy.defaultConfig `jumpPremium;`stepYears;:;1f%252f];
/ Closed: tiny jump intensity -> jumpPx very close to bsPx, threshold huge
closedJump:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(
    0.20;0.01;-0.05;0.05;0.02;0f);
closedCfg:@[baseCfg;(`jumpModelParams;`premiumThreshold);:;(closedJump;5f)];
closedBundle:.strategy.runAndSummarize[`jumpPremium;trade;pathTbl;bsModel;fdmCfg;closedCfg];
.testutil.assertTrue[not closedBundle[`summary;`gateOpen];"gate closed when |premium| < threshold"];
.testutil.assertTrue[0f=closedBundle[`summary;`totalPnl];"totalPnl=0 when flat"];
/ Open: large jump intensity raises premium dramatically; small threshold
openJump:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(
    0.20;3f;-0.05;0.30;0.02;0f);
openCfg:@[baseCfg;(`jumpModelParams;`premiumThreshold);:;(openJump;0.10)];
openBundle:.strategy.runAndSummarize[`jumpPremium;trade;pathTbl;bsModel;fdmCfg;openCfg];
.testutil.assertTrue[openBundle[`summary;`gateOpen];"gate open with large premium and low threshold"];
.testutil.assertTrue[(openBundle[`summary;`tradeSide])=`short;"auto -> short when jump model > BS"];
/ Force-long override
forceLongCfg:@[openCfg;`direction;:;`forceLong];
flBundle:.strategy.runAndSummarize[`jumpPremium;trade;pathTbl;bsModel;fdmCfg;forceLongCfg];
.testutil.assertTrue[(flBundle[`summary;`tradeSide])=`long;"forceLong overrides auto"];
-1 "PASS test_strategy_jump_premium_gate: closedPnl=",string[closedBundle[`summary;`totalPnl]];
