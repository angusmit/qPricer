\l core/init.q
/ hedgeDelta=0b: hedgePosition stays 0, hedgePnl stays 0, hedgeTrade stays 0.
/ __hedgeStep is bypassed but the portfolio-value attribution identity still holds.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;15;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CR_UN;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);

stratCfg:.strategy.defaultConfig `calendarRoll;
stratCfg:@[stratCfg;(`frontTenorYears;`backTenorYears;`rollThresholdYears;`hedgeDelta;`stepYears;`txnCostRate;`financingRate);:;(20f%252f;60f%252f;0.005;0b;1f%252f;0.001;0.02)];

bundle:.strategy.runAndSummarize[`calendarRoll;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;

.testutil.assertTrue[all 0f=resultTbl`hedgePosition;"hedgePosition stays 0 when hedgeDelta=0b"];
.testutil.assertTrue[all 0f=resultTbl`hedgeTrade;"hedgeTrade stays 0"];
.testutil.assertTrue[all 0f=resultTbl`hedgePnl;"hedgePnl stays 0"];

.testutil.assertNear[summary`hedgePnlTotal;0f;1e-12;"hedgePnlTotal = 0"];
.testutil.assertTrue[(summary`financingTotal)<>0f;"financingTotal non-zero (positive financingRate on the cash position)"];

identity:(resultTbl`stepPnl)-(((resultTbl`positionPnl)+(resultTbl`rollPnl)+(resultTbl`hedgePnl)+resultTbl`financingPnl)-resultTbl`txnCost);
maxAbsRes:max abs identity;
.testutil.assertTrue[maxAbsRes<1e-8;"identity holds unhedged within 1e-8"];

-1 "PASS test_strategy_calendar_roll_unhedged: maxIdentityResidual=",string[maxAbsRes],", financingTotal=",string[summary`financingTotal];
