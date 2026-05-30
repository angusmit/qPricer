\l core/init.q
pathCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.05 0.10 0.25 0.50f;
    `simple;
    `spot0`drift`volatility`contango`riskFreeRate!(60f;0f;0.20;1f;0.02);
    6;
    1f%252f;
    23);
bundle:.strategy.path.fromFuturesCurve pathCfg;
pathTbl:bundle`frontPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CC_S;`OIL;`equityOption;`european;`call;60f;0.50;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `commodityCalendar;
stratCfg:@[stratCfg;(`curveBundle;`nearTenor;`farTenor;`rollTriggerTenor;`stepYears);:;(bundle;0.10;0.50;0.25;1f%252f)];
runBundle:.strategy.runAndSummarize[`commodityCalendar;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:runBundle`result;
summary:runBundle`summary;
requiredCols:`stepIndex`stepDate`spot`nearFutures`farFutures`nearEntry`farEntry`positionValue`spreadValue`rollEvents`txnCost`positionPnl`rollPnl`financingPnl`stepPnl`cumulativePnl`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"commodityCalendar schema"];
.testutil.assertTrue[6=count resultTbl;"6 rows"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];
.testutil.assertTrue[(summary`nearTenor)=0.10;"nearTenor reported"];
.testutil.assertTrue[(summary`farTenor)=0.50;"farTenor reported"];
/ Initial position value = 0 (just entered, no markup)
.testutil.assertTrue[1e-12>abs first resultTbl`positionValue;"position value 0 at entry (futures mark-to-market)"];
-1 "PASS test_strategy_commodity_calendar: spreadEntry=",string[first resultTbl`spreadValue],", totalPnl=",string[summary`totalPnl];
