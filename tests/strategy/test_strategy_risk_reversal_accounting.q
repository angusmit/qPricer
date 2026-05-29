\l lib/init.q
/ Independent revaluation: PV computed from state components each step; deltaPV must
/ equal stepPnl from the result table within 1e-8. Catches cash-bookkeeping or
/ leg-mark update bugs that the flow-vs-flow assertion would miss.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;5;1f%252f;0.02;0f;31);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `RR_AC;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `riskReversal;
stratCfg:@[stratCfg;(`skewSlope;`fairSkew;`skewMargin;`stepYears;`txnCostRate;`financingRate);:;(-0.6;-0.2;0.05;1f%252f;0.001;0.02)];

pathRows:0!pathTbl;
firstStep:first pathRows;
remainingRows:1_pathRows;
initialState:.strategy.riskReversal.init[trade;firstStep;bsModel;fdmCfg;stratCfg];
stepFnLocal:.strategy.riskReversal.step[;;trade;bsModel;fdmCfg;stratCfg];
foldedStates:stepFnLocal\[initialState;remainingRows];
allStates:enlist[initialState],foldedStates;
spots:pathRows`spot;

pvFn:{[s;spotVal] .strategy.__portfolioValue[s`cash;s`prevPositionValue;s`hedgePosition;spotVal]};
pvSeries:pvFn'[allStates;spots];
deltaPV:1_(pvSeries-prev pvSeries);

bundle:.strategy.runAndSummarize[`riskReversal;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
stepPnlsAfterInit:1_resultTbl`stepPnl;
independentResidual:max abs deltaPV-stepPnlsAfterInit;
.testutil.assertTrue[independentResidual<1e-8;"deltaPV matches stepPnl within 1e-8 (independent revaluation)"];

flowSumPerRow:((resultTbl`positionPnl)+(resultTbl`hedgePnl)+resultTbl`financingPnl)-resultTbl`txnCost;
flowMaxAbsRes:max abs (resultTbl`stepPnl)-flowSumPerRow;
.testutil.assertTrue[flowMaxAbsRes<1e-8;"flow-bucket identity holds within 1e-8"];

-1 "PASS test_strategy_risk_reversal_accounting: independentResidual=",string[independentResidual],", flowResidual=",string[flowMaxAbsRes];
