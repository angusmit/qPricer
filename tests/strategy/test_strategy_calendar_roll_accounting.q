\l core/init.q
/ Independent-revaluation accounting test. Run calendarRoll via a scan that exposes
/ all per-step states. At each step compute portfolio value INDEPENDENTLY from
/ state components (cash + legMarkSum + hedgePosition*spot) and assert
/ deltaPV == positionPnl + rollPnl + hedgePnl + financingPnl - txnCost within 1e-8.
/ Both sides come from different code paths so the test catches cash-bookkeeping
/ bugs, leg-mark update bugs, or any sign error in the attribution.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.30;25;1f%252f;0.03;0f;19);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CR_AC;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);

stratCfg:.strategy.defaultConfig `calendarRoll;
stratCfg:@[stratCfg;(`frontTenorYears;`backTenorYears;`rollThresholdYears;`stepYears;`txnCostRate;`financingRate);:;(15f%252f;45f%252f;0.005;1f%252f;0.001;0.02)];

pathRows:0!pathTbl;
firstStep:first pathRows;
remainingRows:1_pathRows;
initialState:.strategy.calendarRoll.init[trade;firstStep;bsModel;fdmCfg;stratCfg];
stepFnLocal:.strategy.calendarRoll.step[;;trade;bsModel;fdmCfg;stratCfg];
foldedStates:stepFnLocal\[initialState;remainingRows];
allStates:enlist[initialState],foldedStates;

spots:pathRows`spot;
pvFn:{[s;spotVal] .strategy.__portfolioValue[s`cash;s`prevPositionValue;s`hedgePosition;spotVal]};
pvSeries:pvFn'[allStates;spots];
deltaPV:1_(pvSeries-prev pvSeries);

bundle:.strategy.runAndSummarize[`calendarRoll;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
stepPnlsAfterInit:1_resultTbl`stepPnl;

independentResidual:max abs deltaPV-stepPnlsAfterInit;
.testutil.assertTrue[independentResidual<1e-8;"deltaPV (independent revaluation) matches stepPnl (flow attribution) within 1e-8"];

flowSumPerRow:((resultTbl`positionPnl)+(resultTbl`rollPnl)+(resultTbl`hedgePnl)+resultTbl`financingPnl)-resultTbl`txnCost;
flowMaxAbsRes:max abs (resultTbl`stepPnl)-flowSumPerRow;
.testutil.assertTrue[flowMaxAbsRes<1e-8;"flow-bucket identity holds within 1e-8"];

rollRows:resultTbl where (resultTbl`rollEvents)>0;
nonRollRows:resultTbl where (resultTbl`rollEvents)=0;
.testutil.assertTrue[0<count rollRows;"at least one roll row in this path"];
.testutil.assertTrue[all 0f=nonRollRows`rollPnl;"rollPnl is exactly zero on non-roll steps"];

summary:bundle`summary;
.testutil.assertTrue[`OK=summary`status;"summary status OK"];
.testutil.assertNear[summary`totalPnl;sum resultTbl`stepPnl;1e-10;"summary totalPnl = sum stepPnl"];

-1 "PASS test_strategy_calendar_roll_accounting: independentResidual=",string[independentResidual],", flowResidual=",string[flowMaxAbsRes],", rollRows=",string[count rollRows];
