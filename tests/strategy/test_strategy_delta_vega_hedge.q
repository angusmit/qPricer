\l lib/init.q
/ deltaVegaHedge structural + neutrality test. After each rebalance, the residual
/ vega (bookVega + vegaHedgeUnits * hedgeOptionVega) should be at machine epsilon.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;7;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `DVH_S;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `deltaVegaHedge;
stratCfg:@[stratCfg;`stepYears;:;1f%252f];

bundle:.strategy.runAndSummarize[`deltaVegaHedge;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;

requiredCols:`stepIndex`stepDate`spot`volatility`bookPrice`hedgeOptionPrice`vegaHedgeUnits`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`bookPnl`vegaHedgePnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`residualVega`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"deltaVegaHedge schema"];
.testutil.assertTrue[7=count resultTbl;"7 rows"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];
.testutil.assertTrue[all 1e-10>abs resultTbl`residualVega;"residual vega within machine epsilon every step (vega-neutral)"];
.testutil.assertTrue[(summary`maxAbsResidualVega)<1e-10;"summary maxAbsResidualVega very small"];

.testutil.assertTrue[(summary`numVegaRehedges)>=1;"vega rebalances counted"];
.testutil.assertTrue[`OK=summary`status;"summary OK"];

bookPnlAbs:abs summary`bookPnlTotal;
totalPnlAbs:abs summary`totalPnl;
.testutil.assertTrue[totalPnlAbs<0.5*bookPnlAbs;"delta+vega neutral book has totalPnl much smaller than book pnl alone"];

-1 "PASS test_strategy_delta_vega_hedge: maxAbsResidualVega=",string[summary`maxAbsResidualVega],", totalPnl=",string[summary`totalPnl],", bookPnlTotal=",string[summary`bookPnlTotal];
