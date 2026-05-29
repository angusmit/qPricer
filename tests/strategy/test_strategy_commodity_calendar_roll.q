\l lib/init.q
/ Use a long-enough path so the near tenor (0.05 years ~ 12.6 trading days) expires.
pathCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.05 0.10 0.25 0.50f;
    `simple;
    `spot0`drift`volatility`contango`riskFreeRate!(60f;0f;0.20;1f;0.02);
    20;
    1f%252f;
    23);
bundle:.strategy.path.fromFuturesCurve pathCfg;
pathTbl:bundle`frontPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CC_R;`OIL;`equityOption;`european;`call;60f;0.50;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `commodityCalendar;
/ nearTenor 0.05 yrs ~ step 13 will trigger roll; rollTriggerTenor 0.10 used as next near
stratCfg:@[stratCfg;(`curveBundle;`nearTenor;`farTenor;`rollTriggerTenor;`stepYears);:;(bundle;0.05;0.50;0.10;1f%252f)];
runBundle:.strategy.runAndSummarize[`commodityCalendar;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:runBundle`result;
summary:runBundle`summary;
.testutil.assertTrue[(summary`totalRolls)>0;"at least one roll event"];
rollRow:resultTbl where (resultTbl`rollEvents)>0;
.testutil.assertTrue[0<count rollRow;"a row marks the roll"];
.testutil.assertTrue[(first rollRow`nearEntry)<>first resultTbl`nearEntry;"nearEntry resets after roll"];
-1 "PASS test_strategy_commodity_calendar_roll: totalRolls=",string[summary`totalRolls];
