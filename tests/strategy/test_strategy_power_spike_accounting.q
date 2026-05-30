\l core/init.q
pathCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.05 0.1 0.25 0.5f;
    `mrjump;
    `initialLogPrice`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility`contango`riskFreeRate!(
        log 50f;2f;log 50f;0.30;30f;0.10;0.20;0f;0.02);
    5;
    1f%252f;
    11);
bundle:.strategy.path.fromFuturesCurve pathCfg;
pathTbl:bundle`frontPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `PSC_AC;`PWR;`equityOption;`european;`call;50f;0.10;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `powerSpikeCapture;
stratCfg:@[stratCfg;(`curveBundle;`stepYears;`financingRate;`txnCostRate);:;(bundle;1f%252f;0.02;0.001)];
pathRows:0!pathTbl;
firstStep:first pathRows;
remainingRows:1_pathRows;
initialState:.strategy.powerSpikeCapture.init[trade;firstStep;bsModel;fdmCfg;stratCfg];
stepFnLocal:.strategy.powerSpikeCapture.step[;;trade;bsModel;fdmCfg;stratCfg];
foldedStates:stepFnLocal\[initialState;remainingRows];
allStates:enlist[initialState],foldedStates;
spots:pathRows`spot;
pvFn:{[s;sp] .strategy.__portfolioValue[s`cash;s`prevPositionValue;s`hedgePosition;sp]};
pvSeries:pvFn'[allStates;spots];
deltaPV:1_(pvSeries-prev pvSeries);
bundleRun:.strategy.runAndSummarize[`powerSpikeCapture;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundleRun`result;
stepPnls:1_resultTbl`stepPnl;
independentResidual:max abs deltaPV-stepPnls;
.testutil.assertTrue[independentResidual<1e-8;"powerSpike deltaPV matches stepPnl within 1e-8"];
-1 "PASS test_strategy_power_spike_accounting: indepResidual=",string[independentResidual];
