\l core/init.q
pathCfg:`names`weights`spot0`vols`correlationMatrix`steps`stepYears`riskFreeRate`dividendYields`seed!(
    `A`B;0.5 0.5;100 100f;0.20 0.30;(1 0.3f;0.3 1f);4;1f%252f;0.02;0 0f;13);
bundle:.strategy.path.fromCorrelated pathCfg;
pathTbl:bundle`indexPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `DS_AC;`IDX;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `dispersion;
stratCfg:@[stratCfg;(`pathBundle;`indexVol;`stepYears;`financingRate;`txnCostRate);:;(bundle;0.30;1f%252f;0.02;0.001)];
pathRows:0!pathTbl;
firstStep:first pathRows;
remainingRows:1_pathRows;
initialState:.strategy.dispersion.init[trade;firstStep;bsModel;fdmCfg;stratCfg];
stepFnLocal:.strategy.dispersion.step[;;trade;bsModel;fdmCfg;stratCfg];
foldedStates:stepFnLocal\[initialState;remainingRows];
allStates:enlist[initialState],foldedStates;
namesList:bundle`names;
pathTables:bundle`pathTables;
nNames:count namesList;
/ Independent PV at step k = state.cash + state.prevPositionValue + sum(hedgePos_i * spot_i_at_step_k)
pvFn:{[s]
    hedgePositions:exec hedgePosition from s`hedgeBook;
    spotsAtStep:s`prevSpots;
    .strategy.__portfolioValueMulti[s`cash;s`prevPositionValue;hedgePositions;spotsAtStep]};
pvSeries:pvFn each allStates;
deltaPV:1_(pvSeries-prev pvSeries);
bundleRun:.strategy.runAndSummarize[`dispersion;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundleRun`result;
stepPnls:1_resultTbl`stepPnl;
independentResidual:max abs deltaPV-stepPnls;
.testutil.assertTrue[independentResidual<1e-8;"dispersion deltaPV matches stepPnl within 1e-8"];
-1 "PASS test_strategy_dispersion_accounting: indepResidual=",string[independentResidual];
