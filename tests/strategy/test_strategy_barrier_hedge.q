\l lib/init.q
/ A path that stays comfortably below the up-and-out barrier - option survives.
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.15;6;1f%252f;0.02;0f;3);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel!(
    `BH_S;`X;`equityOption;`european;`call;100f;0.25;1f;`upAndOut;150f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;30;500;0f;300f;`linear;1b;0b);
stratCfg:.strategy.defaultConfig `barrierHedge;
stratCfg:@[stratCfg;`stepYears;:;1f%252f];
bundle:.strategy.runAndSummarize[`barrierHedge;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;
requiredCols:`stepIndex`stepDate`spot`barrierLevel`distanceToBarrier`knockedOut`optionPrice`delta`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"barrierHedge schema"];
.testutil.assertTrue[6=count resultTbl;"6 rows"];
.testutil.assertTrue[not summary`knockedOut;"barrier not breached on benign path"];
.testutil.assertTrue[all not resultTbl`knockedOut;"all steps live"];
.testutil.assertTrue[all 0f<resultTbl`optionPrice;"option prices positive on live steps"];
.testutil.assertTrue[all (resultTbl`barrierLevel)=150f;"barrierLevel reported"];
.testutil.assertTrue[all (resultTbl`distanceToBarrier)>0f;"distance to barrier positive on live path"];
-1 "PASS test_strategy_barrier_hedge: maxAbsHedgeTrade=",string[summary`maxAbsHedgeTrade],", totalPnl=",string[summary`totalPnl];
