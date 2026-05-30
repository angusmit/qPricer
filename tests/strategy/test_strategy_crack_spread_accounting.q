\l core/init.q
/ Independent-revaluation accounting check for crackSpread (PV = cash + optionMark).
pathCfg:`names`correlationMatrix`spot0s`drifts`vols`contangos`tenors`steps`stepYears`riskFreeRate`seed!(
    `product`crude;
    (1 0.6f;0.6 1f);
    95 80f;
    0 0f;
    0.35 0.30f;
    1 0.8f;
    0.08 0.25 0.5f;
    6;
    1f%252f;
    0.03;
    88);
bundle:.strategy.path.fromCorrelatedCurves pathCfg;
pathTbl:(bundle`curves)[`product]`frontPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CRK_AC;`PRODUCT;`equityOption;`european;`call;0f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `crackSpread;
stratCfg:@[stratCfg;(`curveBundle;`crackRatio;`vol1;`vol2;`correlation;`expiry;`stepYears;`financingRate;`txnCostRate);:;(bundle;1f;0.35;0.30;0.6;0.25;1f%252f;0.03;0.0005)];

pathRows:0!pathTbl;
firstStep:first pathRows;
remainingRows:1_pathRows;
normSpec:.strategy.crackSpread.__normSpec stratCfg;
initialState:.strategy.spreadOption.coreInit[trade;firstStep;stratCfg;normSpec];
stepFnLocal:.strategy.spreadOption.coreStep[;;stratCfg];
foldedStates:stepFnLocal\[initialState;remainingRows];
allStates:enlist[initialState],foldedStates;

pvFn:{[s] .strategy.__portfolioValue[s`cash;s`optionMark;0f;0f]};
pvSeries:pvFn each allStates;
deltaPV:1_(pvSeries-prev pvSeries);

bundleRun:.strategy.runAndSummarize[`crackSpread;trade;pathTbl;bsModel;fdmCfg;stratCfg];
stepPnls:1_(bundleRun`result)`stepPnl;
independentResidual:max abs deltaPV-stepPnls;
.testutil.assertTrue[independentResidual<1e-8;"crackSpread deltaPV matches stepPnl within 1e-8"];

-1 "PASS test_strategy_crack_spread_accounting: indepResidual=",string independentResidual;
