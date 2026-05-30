\l core/init.q
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.15;7;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `LV_S;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `longVol;
stratCfg:@[stratCfg;(`forecastVol;`entryMargin;`stepYears);:;(0.30;0.02;1f%252f)];
bundle:.strategy.runAndSummarize[`longVol;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;
requiredCols:`stepIndex`stepDate`spot`volatility`callPrice`putPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`premiumPaid`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"longVol schema"];
.testutil.assertTrue[summary`gateOpen;"gate open: forecast 0.30 > implied 0.15 + margin 0.02"];
.testutil.assertTrue[(summary`premiumPaid)>0f;"premium paid > 0"];
.testutil.assertTrue[(summary`theoreticalGammaPnlTotal)>0f;"long gamma -> positive theoretical gamma PnL in moving market"];
.testutil.assertTrue[(summary`thetaPnlTotal)<0f;"long theta -> negative time decay"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];
-1 "PASS test_strategy_long_vol: premium=",string[summary`premiumPaid],", gammaTot=",string[summary`theoreticalGammaPnlTotal];
