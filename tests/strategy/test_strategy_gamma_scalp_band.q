\l core/init.q
/ Band-mode rebalancing rebalances less often than interval-mode on the same path,
/ and with a positive transaction-cost rate that translates into lower txnCostTotal.

pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.40;15;1f%252f;0.05;0f;13);
pathTbl:.strategy.path.fromSynthetic pathCfg;

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `T_BAND;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;150;200;0f;300f;`linear;1b;1b);

baseCfg:.strategy.defaultConfig `gammaScalp;
baseCfg:@[baseCfg;(`stepYears;`txnCostRate);:;(1f%252f;0.001)];

intervalCfg:@[baseCfg;(`rebalanceMode;`rebalanceInterval);:;(`interval;1)];
bandCfg:@[baseCfg;(`rebalanceMode;`deltaBand);:;(`band;0.10)];

intervalBundle:.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;intervalCfg];
bandBundle:.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;bandCfg];

intervalRebals:intervalBundle[`summary;`numRebalances];
bandRebals:bandBundle[`summary;`numRebalances];
intervalCost:intervalBundle[`summary;`txnCostTotal];
bandCost:bandBundle[`summary;`txnCostTotal];

.testutil.assertTrue[bandRebals<=intervalRebals;"band mode rebalances <= interval mode rebalances"];
.testutil.assertTrue[bandCost<=intervalCost;"band mode txnCostTotal <= interval mode txnCostTotal"];

.testutil.assertTrue[(intervalCost)>0f;"interval mode incurs transaction cost"];
.testutil.assertTrue[(intervalBundle`summary)[`status]=`OK;"interval bundle status OK"];
.testutil.assertTrue[(bandBundle`summary)[`status]=`OK;"band bundle status OK"];

-1 "PASS test_strategy_gamma_scalp_band: intervalRebals=",string[intervalRebals],", bandRebals=",string[bandRebals],", intervalCost=",string[intervalCost],", bandCost=",string[bandCost];
