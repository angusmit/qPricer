\l core/init.q
/ Independent-revaluation accounting check for sparkSpread. Reconstruct the
/ portfolio value PV = cash + optionMark (futures hedges carry zero mark; their
/ P&L is realized to cash) from state components and assert deltaPV == stepPnl.
pathCfg:`names`correlationMatrix`spot0s`drifts`vols`contangos`tenors`steps`stepYears`riskFreeRate`seed!(
    `power`gas;
    (1 0.4f;0.4 1f);
    50 3.5f;
    0 0f;
    0.45 0.35f;
    0.5 0.1f;
    0.08 0.25 0.5f;
    6;
    1f%252f;
    0.03;
    77);
bundle:.strategy.path.fromCorrelatedCurves pathCfg;
pathTbl:(bundle`curves)[`power]`frontPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `SPK_AC;`POWER;`equityOption;`european;`call;0f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `sparkSpread;
stratCfg:@[stratCfg;(`curveBundle;`heatRate;`vol1;`vol2;`correlation;`expiry;`stepYears;`financingRate;`txnCostRate);:;(bundle;8f;0.45;0.35;0.4;0.25;1f%252f;0.03;0.0005)];

pathRows:0!pathTbl;
firstStep:first pathRows;
remainingRows:1_pathRows;
normSpec:.strategy.sparkSpread.__normSpec stratCfg;
initialState:.strategy.spreadOption.coreInit[trade;firstStep;stratCfg;normSpec];
stepFnLocal:.strategy.spreadOption.coreStep[;;stratCfg];
foldedStates:stepFnLocal\[initialState;remainingRows];
allStates:enlist[initialState],foldedStates;

/ Independent PV: cash + optionMark (single-asset helper with zero hedge term).
pvFn:{[s] .strategy.__portfolioValue[s`cash;s`optionMark;0f;0f]};
pvSeries:pvFn each allStates;
deltaPV:1_(pvSeries-prev pvSeries);

bundleRun:.strategy.runAndSummarize[`sparkSpread;trade;pathTbl;bsModel;fdmCfg;stratCfg];
stepPnls:1_(bundleRun`result)`stepPnl;
independentResidual:max abs deltaPV-stepPnls;
.testutil.assertTrue[independentResidual<1e-8;"sparkSpread deltaPV matches stepPnl within 1e-8"];

-1 "PASS test_strategy_spark_spread_accounting: indepResidual=",string independentResidual;
