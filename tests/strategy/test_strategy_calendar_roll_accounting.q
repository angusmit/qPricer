\l lib/init.q
/ Portfolio-value attribution identity: stepPnl == positionPnl + rollPnl + hedgePnl
/ + financingPnl - txnCost for every row, within 1e-8. Also: rollPnl is non-zero on
/ rows with rollEvents>0 (typically) and zero elsewhere; positionPnl is the surviving-
/ leg mark change.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.30;25;1f%252f;0.03;0f;19);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CR_AC;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);

stratCfg:.strategy.defaultConfig `calendarRoll;
stratCfg:@[stratCfg;(`frontTenorYears;`backTenorYears;`rollThresholdYears;`stepYears;`txnCostRate;`financingRate);:;(15f%252f;45f%252f;0.005;1f%252f;0.001;0.02)];

bundle:.strategy.runAndSummarize[`calendarRoll;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;

identity:(resultTbl`stepPnl)-(((resultTbl`positionPnl)+(resultTbl`rollPnl)+(resultTbl`hedgePnl)+resultTbl`financingPnl)-resultTbl`txnCost);
maxAbsRes:max abs identity;
.testutil.assertTrue[maxAbsRes<1e-8;"identity stepPnl == positionPnl + rollPnl + hedgePnl + financingPnl - txnCost (max residual < 1e-8)"];

rollRows:resultTbl where (resultTbl`rollEvents)>0;
nonRollRows:resultTbl where (resultTbl`rollEvents)=0;
.testutil.assertTrue[0<count rollRows;"at least one roll row in this path"];
.testutil.assertTrue[0<count nonRollRows;"at least one non-roll row"];
.testutil.assertTrue[all 0f=nonRollRows`rollPnl;"rollPnl is exactly zero on non-roll steps"];

summary:bundle`summary;
.testutil.assertTrue[`OK=summary`status;"summary status OK"];
.testutil.assertNear[summary`totalPnl;sum resultTbl`stepPnl;1e-10;"summary totalPnl = sum stepPnl"];

-1 "PASS test_strategy_calendar_roll_accounting: maxIdentityResidual=",string[maxAbsRes],", rollRows=",string[count rollRows];
