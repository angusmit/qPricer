\l core/init.q
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;7;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `IC_S;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `ironCondor;
stratCfg:@[stratCfg;`stepYears;:;1f%252f];
bundle:.strategy.runAndSummarize[`ironCondor;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;
requiredCols:`stepIndex`stepDate`spot`shortPutPrice`longPutPrice`shortCallPrice`longCallPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"ironCondor schema"];
.testutil.assertTrue[7=count resultTbl;"7 rows"];
.testutil.assertTrue[summary`gateOpen;"gate open"];
.testutil.assertTrue[(summary`netCredit)>0f;"net credit positive"];
.testutil.assertTrue[all 0f<resultTbl`shortPutPrice;"short put prices positive"];
.testutil.assertTrue[all 0f<resultTbl`longPutPrice;"long put prices positive"];
.testutil.assertTrue[all 0f<resultTbl`shortCallPrice;"short call prices positive"];
.testutil.assertTrue[all 0f<resultTbl`longCallPrice;"long call prices positive"];
.testutil.assertTrue[(summary`maxProfit)=summary`netCredit;"maxProfit = netCredit"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];
-1 "PASS test_strategy_iron_condor: netCredit=",string[summary`netCredit],", maxLoss=",string[summary`maxLoss];
