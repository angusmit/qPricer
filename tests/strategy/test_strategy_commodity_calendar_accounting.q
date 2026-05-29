\l lib/init.q
pathCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.05 0.10 0.25 0.50f;
    `simple;
    `spot0`drift`volatility`contango`riskFreeRate!(60f;0f;0.20;1f;0.02);
    8;
    1f%252f;
    23);
bundle:.strategy.path.fromFuturesCurve pathCfg;
pathTbl:bundle`frontPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CC_AC;`OIL;`equityOption;`european;`call;60f;0.50;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `commodityCalendar;
stratCfg:@[stratCfg;(`curveBundle;`nearTenor;`farTenor;`rollTriggerTenor;`stepYears;`financingRate;`txnCostRate);:;(bundle;0.10;0.50;0.25;1f%252f;0.02;0.001)];
pathRows:0!pathTbl;
firstStep:first pathRows;
remainingRows:1_pathRows;
initialState:.strategy.commodityCalendar.init[trade;firstStep;bsModel;fdmCfg;stratCfg];
stepFnLocal:.strategy.commodityCalendar.step[;;trade;bsModel;fdmCfg;stratCfg];
foldedStates:stepFnLocal\[initialState;remainingRows];
allStates:enlist[initialState],foldedStates;
spots:pathRows`spot;
pvFn:{[s;sp] .strategy.__portfolioValue[s`cash;s`prevPositionValue;s`hedgePosition;sp]};
pvSeries:pvFn'[allStates;spots];
deltaPV:1_(pvSeries-prev pvSeries);
bundleRun:.strategy.runAndSummarize[`commodityCalendar;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundleRun`result;
stepPnls:1_resultTbl`stepPnl;
independentResidual:max abs deltaPV-stepPnls;
.testutil.assertTrue[independentResidual<1e-8;"commodityCalendar deltaPV matches stepPnl within 1e-8"];
-1 "PASS test_strategy_commodity_calendar_accounting: indepResidual=",string[independentResidual];
