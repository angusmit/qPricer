\l core/init.q
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;7;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CT_S;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `collarTailHedge;
stratCfg:@[stratCfg;`stepYears;:;1f%252f];
bundle:.strategy.runAndSummarize[`collarTailHedge;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;
requiredCols:`stepIndex`stepDate`spot`mode`underlyingValue`putValue`callValue`positionValue`netDelta`txnCost`positionPnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`putUnits`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"collarTailHedge schema"];
.testutil.assertTrue[7=count resultTbl;"7 rows"];
.testutil.assertTrue[`collar=summary`mode;"mode = collar"];
.testutil.assertTrue[all 0f<resultTbl`underlyingValue;"underlying value positive"];
.testutil.assertTrue[all 0f<resultTbl`putValue;"put value positive"];
.testutil.assertTrue[(summary`protectionFloor)<100f;"protectionFloor below spot0"];
.testutil.assertTrue[(summary`cap)>100f;"cap above spot0"];
-1 "PASS test_strategy_collar_tail_hedge: protectionFloor=",string[summary`protectionFloor],", cap=",string[summary`cap];
