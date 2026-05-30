\l core/init.q
/ Delta + vega neutrality assertion. After each rebalance, netDelta and residualVega
/ should both be near zero in absolute value.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.30;9;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `DVH_N;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `deltaVegaHedge;
stratCfg:@[stratCfg;(`stepYears;`rebalanceInterval);:;(1f%252f;1)];

bundle:.strategy.runAndSummarize[`deltaVegaHedge;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
postRebalanceRows:1_resultTbl;

.testutil.assertTrue[all 1e-10>abs postRebalanceRows`residualVega;"residual vega ~ 0 after each vega rehedge"];

postHedgeResidualDelta:(postRebalanceRows`netDelta)+postRebalanceRows`hedgePosition;
.testutil.assertTrue[all 1e-10>abs postHedgeResidualDelta;"residual delta (netDelta + hedgePosition) ~ 0 after each delta hedge"];

initRowSelect:resultTbl 0;
.testutil.assertTrue[1e-10>abs initRowSelect`residualVega;"residual vega ~ 0 at init too"];

-1 "PASS test_strategy_delta_vega_hedge_neutrality: maxResVega=",string[max abs resultTbl`residualVega],", maxNetDelta=",string[max abs resultTbl`netDelta];
