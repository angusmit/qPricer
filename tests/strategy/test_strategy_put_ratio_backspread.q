\l lib/init.q
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.30;7;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `PR_S;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:@[.strategy.defaultConfig `putRatioBackspread;`stepYears;:;1f%252f];
bundle:.strategy.runAndSummarize[`putRatioBackspread;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;
requiredCols:`stepIndex`stepDate`spot`shortStrike`longStrike`shortPutPrice`longPutPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"putRatioBackspread schema"];
.testutil.assertTrue[7=count resultTbl;"7 rows"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];
.testutil.assertTrue[(summary`shortStrike)>summary`longStrike;"shortStrike > longStrike"];
.testutil.assertTrue[2f=summary`ratioN;"ratioN = 2 by default"];
.testutil.assertTrue[(summary`theoreticalGammaPnlTotal)>0f;"net long ratio puts -> positive gamma in moving market"];
-1 "PASS test_strategy_put_ratio_backspread: shortStrike=",string[summary`shortStrike],", longStrike=",string[summary`longStrike];
