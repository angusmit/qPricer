\l lib/init.q
/ A path long enough to expire the front leg at least once triggers totalRolls >= 1,
/ activeLegCount returns to 2 after a roll (rollBackLeg=1b keeps constant maturity, but
/ here we test the typical case where front expires while back is still alive), and the
/ front's remaining time resets upward at the roll step.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;25;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CR_EV;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);

stratCfg:.strategy.defaultConfig `calendarRoll;
stratCfg:@[stratCfg;(`frontTenorYears;`backTenorYears;`rollThresholdYears;`rollBackLeg;`stepYears);:;(15f%252f;45f%252f;0.005;1b;1f%252f)];

bundle:.strategy.runAndSummarize[`calendarRoll;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;

.testutil.assertTrue[(summary`totalRolls)>=1;"at least one roll event in 25-step path with 15-day front"];
.testutil.assertTrue[`OK=summary`status;"summary status OK"];

rollSteps:resultTbl where (resultTbl`rollEvents)>0;
.testutil.assertTrue[0<count rollSteps;"at least one row marks rollEvents>0"];

firstRollRow:first rollSteps;
.testutil.assertTrue[(firstRollRow`activeLegCount)>=2;"activeLegCount >= 2 after a roll under rollBackLeg=1b"];

priorIdx:firstRollRow[`stepIndex]-1;
priorRow:resultTbl first where (resultTbl`stepIndex)=priorIdx;
.testutil.assertTrue[(firstRollRow`frontRemainingTime)>priorRow`frontRemainingTime;"frontRemainingTime resets upward after roll"];

postRollRow:resultTbl first where (resultTbl`stepIndex)=1+firstRollRow`stepIndex;
.testutil.assertTrue[(postRollRow`frontRemainingTime)<firstRollRow`frontRemainingTime;"front decays again after roll"];

-1 "PASS test_strategy_calendar_roll_event: totalRolls=",string[summary`totalRolls],", firstRollStep=",string[firstRollRow`stepIndex];
