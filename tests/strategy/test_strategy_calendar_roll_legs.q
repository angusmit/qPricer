\l core/init.q
/ calendarRoll.init builds exactly two legs for a longCalendar: front (short, shorter
/ expiry) + back (long, longer expiry). For shortCalendar the sides are flipped.
/ The legs table follows the documented schema; both legs share the input strike and
/ option type from config.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.20;3;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CR_LEG;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);

stratLong:.strategy.defaultConfig `calendarRoll;
stratLong:@[stratLong;(`frontTenorYears;`backTenorYears;`rollThresholdYears;`stepYears);:;(0.10;0.30;0.005;1f%252f)];

stateLong:.strategy.calendarRoll.init[trade;first 0!pathTbl;bsModel;fdmCfg;stratLong];
legs:stateLong`legs;

requiredCols:`legId`legRole`optionType`strike`side`units`expiryTimeYears`entryTimeYears`entrySpot`entryPrice`currentMark`currentDelta`currentGamma`currentTheta;
.testutil.assertTableColumns[legs;requiredCols;"legs table schema"];
.testutil.assertTrue[2=count legs;"init builds exactly 2 legs"];
.testutil.assertTrue[(`back`front)~asc legs`legRole;"roles are front and back"];
.testutil.assertTrue[all (legs`strike)=100f;"both legs share strike from trade"];
.testutil.assertTrue[all (legs`optionType)=`call;"both legs share optionType from config"];

frontLeg:legs first where (legs`legRole)=`front;
backLeg:legs first where (legs`legRole)=`back;
.testutil.assertNear[frontLeg`side;-1f;1e-12;"longCalendar -> front is short"];
.testutil.assertNear[backLeg`side;1f;1e-12;"longCalendar -> back is long"];
.testutil.assertNear[frontLeg`expiryTimeYears;0.10;1e-12;"front expiry = frontTenorYears"];
.testutil.assertNear[backLeg`expiryTimeYears;0.30;1e-12;"back expiry = backTenorYears"];
.testutil.assertTrue[(frontLeg`expiryTimeYears)<backLeg`expiryTimeYears;"front expires before back"];
.testutil.assertTrue[(frontLeg`entryPrice)>0f;"front entry price positive"];
.testutil.assertTrue[(backLeg`entryPrice)>0f;"back entry price positive"];

stratShort:@[stratLong;`spreadType;:;`shortCalendar];
stateShort:.strategy.calendarRoll.init[trade;first 0!pathTbl;bsModel;fdmCfg;stratShort];
legsShort:stateShort`legs;
frontShort:legsShort first where (legsShort`legRole)=`front;
backShort:legsShort first where (legsShort`legRole)=`back;
.testutil.assertNear[frontShort`side;1f;1e-12;"shortCalendar -> front is long"];
.testutil.assertNear[backShort`side;-1f;1e-12;"shortCalendar -> back is short"];

rowEmit:stateLong`rowEmit;
.testutil.assertTrue[`OK=rowEmit`status;"init rowEmit status OK"];
.testutil.assertTrue[2=rowEmit`activeLegCount;"init rowEmit activeLegCount = 2"];
.testutil.assertNear[rowEmit`frontRemainingTime;0.10;1e-12;"init rowEmit frontRemainingTime = front tenor"];

badCfg:@[stratLong;`backTenorYears;:;0.05];
badRes:.[.strategy.calendarRoll.init;(trade;first 0!pathTbl;bsModel;fdmCfg;badCfg);{`ERROR}];
.testutil.assertTrue[badRes~`ERROR;"backTenorYears <= frontTenorYears rejected"];

-1 "PASS test_strategy_calendar_roll_legs: longFrontSide=",string[frontLeg`side],", longBackSide=",string[backLeg`side];
