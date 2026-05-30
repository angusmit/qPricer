\l core/init.q
/ Independent revaluation across a path that includes a knock-out event.
nSteps:5;
spots:100 104 108 112 116f;
pathTbl:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status!(
    til nSteps;
    2024.01.01+til nSteps;
    spots;
    nSteps#0.20;
    nSteps#0.02;
    nSteps#0f;
    nSteps#0Nf;
    nSteps#`OK);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel!(
    `BH_AC;`X;`equityOption;`european;`call;100f;0.25;1f;`upAndOut;110f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;30;500;0f;300f;`linear;1b;0b);
stratCfg:@[.strategy.defaultConfig `barrierHedge;(`stepYears;`financingRate;`txnCostRate);:;(1f%252f;0.02;0.001)];
pathRows:0!pathTbl;
firstStep:first pathRows;
remainingRows:1_pathRows;
initialState:.strategy.barrierHedge.init[trade;firstStep;bsModel;fdmCfg;stratCfg];
stepFnLocal:.strategy.barrierHedge.step[;;trade;bsModel;fdmCfg;stratCfg];
foldedStates:stepFnLocal\[initialState;remainingRows];
allStates:enlist[initialState],foldedStates;
spots:pathRows`spot;
pvFn:{[s;sp] .strategy.__portfolioValue[s`cash;s`prevPositionValue;s`hedgePosition;sp]};
pvSeries:pvFn'[allStates;spots];
deltaPV:1_(pvSeries-prev pvSeries);
bundle:.strategy.runAndSummarize[`barrierHedge;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
stepPnls:1_resultTbl`stepPnl;
independentResidual:max abs deltaPV-stepPnls;
.testutil.assertTrue[independentResidual<1e-8;"deltaPV matches stepPnl within 1e-8 across knock event"];
-1 "PASS test_strategy_barrier_hedge_accounting: indepResidual=",string[independentResidual],", knocked=",string[bundle[`summary;`knockedOut]];
