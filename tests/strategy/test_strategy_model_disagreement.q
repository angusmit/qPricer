\l lib/init.q
/ modelDisagreement structural test.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;7;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `MD;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `modelDisagreement;
stratCfg:@[stratCfg;(`modelAVolBump;`modelBVolBump;`disagreementThreshold;`stepYears);:;(0f;0.10;0.10;1f%252f)];

bundle:.strategy.runAndSummarize[`modelDisagreement;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;

requiredCols:`stepIndex`stepDate`spot`volatility`optionPrice`delta`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`gateOpen`tradeSide`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"modelDisagreement schema"];
.testutil.assertTrue[7=count resultTbl;"7 rows"];
.testutil.assertTrue[summary`gateOpen;"vol bump 0.10 produces disagreement > 0.10"];
.testutil.assertTrue[1=summary`tradeSide;"modelB > modelA -> long"];
.testutil.assertTrue[`OK=summary`status;"summary OK"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all rows OK"];

-1 "PASS test_strategy_model_disagreement: tradeSide=",string[summary`tradeSide],", totalPnl=",string[summary`totalPnl];
