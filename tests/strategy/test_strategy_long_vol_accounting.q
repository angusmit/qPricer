\l core/init.q
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.15;5;1f%252f;0.02;0f;31);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `LV_AC;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
stratCfg:@[.strategy.defaultConfig `longVol;(`forecastVol;`entryMargin;`stepYears;`txnCostRate;`financingRate);:;(0.30;0.02;1f%252f;0.001;0.02)];
pathRows:0!pathTbl;
firstStep:first pathRows;
remainingRows:1_pathRows;
initialState:.strategy.longVol.init[trade;firstStep;bsModel;fdmCfg;stratCfg];
stepFnLocal:.strategy.longVol.step[;;trade;bsModel;fdmCfg;stratCfg];
foldedStates:stepFnLocal\[initialState;remainingRows];
allStates:enlist[initialState],foldedStates;
spots:pathRows`spot;
pvFn:{[s;sp] .strategy.__portfolioValue[s`cash;s`prevPositionValue;s`hedgePosition;sp]};
pvSeries:pvFn'[allStates;spots];
deltaPV:1_(pvSeries-prev pvSeries);
bundle:.strategy.runAndSummarize[`longVol;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
stepPnls:1_resultTbl`stepPnl;
independentResidual:max abs deltaPV-stepPnls;
.testutil.assertTrue[independentResidual<1e-8;"deltaPV matches stepPnl within 1e-8"];
flowSum:((resultTbl`positionPnl)+(resultTbl`hedgePnl)+resultTbl`financingPnl)-resultTbl`txnCost;
flowRes:max abs (resultTbl`stepPnl)-flowSum;
.testutil.assertTrue[flowRes<1e-8;"flow-bucket identity"];
-1 "PASS test_strategy_long_vol_accounting: indepResidual=",string[independentResidual];
