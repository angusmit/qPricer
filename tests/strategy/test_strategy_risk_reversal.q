\l lib/init.q
/ riskReversal structural test: builds 2 wings at +/- wingOffsetPct, per-strike skew
/ vol formula puts the lower vol on the call wing for a negative skew (longCallWing).

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;7;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `RR_S;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `riskReversal;
stratCfg:@[stratCfg;(`skewSlope;`fairSkew;`skewMargin;`stepYears;`wingOffsetPct);:;(-0.6;-0.2;0.05;1f%252f;0.05)];

bundle:.strategy.runAndSummarize[`riskReversal;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;

requiredCols:`stepIndex`stepDate`spot`callStrike`putStrike`callVol`putVol`callPrice`putPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`gateOpen`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"riskReversal result schema"];
.testutil.assertTrue[7=count resultTbl;"7 result rows"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];

firstRow:resultTbl 0;
.testutil.assertNear[firstRow`callStrike;105f;1e-12;"callStrike = spot*1.05"];
.testutil.assertNear[firstRow`putStrike;95f;1e-12;"putStrike = spot*0.95"];
.testutil.assertTrue[firstRow[`callVol]<firstRow`putVol;"negative skew: callVol < putVol"];
.testutil.assertTrue[(firstRow`callVol)>0f;"callVol positive"];
.testutil.assertTrue[(firstRow`putVol)>0f;"putVol positive"];

.testutil.assertTrue[summary`gateOpen;"gate open under these params"];
.testutil.assertTrue[`OK=summary`status;"summary OK"];

-1 "PASS test_strategy_risk_reversal: callStrike=",string[firstRow`callStrike],", callVol=",string[firstRow`callVol],", putVol=",string[firstRow`putVol];
