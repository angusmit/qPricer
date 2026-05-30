\l core/init.q
/ Run delta-hedged gamma scalp on a synthetic path. Structural checks:
/ result is a finite-pnl table with one row per step, summary has the expected
/ totals, and for a long-gamma position in a moving market theoreticalGammaPnlTotal
/ is positive.

pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.30;11;1f%252f;0.05;0f;42);
pathTbl:.strategy.path.fromSynthetic pathCfg;

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `T_GS;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;200;0f;300f;`linear;1b;1b);

stratCfg:.strategy.defaultConfig `gammaScalp;
stratCfg:@[stratCfg;(`rebalanceMode;`rebalanceInterval;`stepYears);:;(`interval;1;1f%252f)];

bundle:.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;

requiredCols:`stepIndex`stepDate`spot`volatility`optionPrice`delta`hedgePosition`hedgeTrade`txnCost`optionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"result schema"];
.testutil.assertTrue[11=count resultTbl;"result has 11 rows (one per path step)"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];
.testutil.assertTrue[all not null resultTbl`cumulativePnl;"cumulativePnl finite"];
.testutil.assertTrue[all 0f<resultTbl`optionPrice;"optionPrice positive"];

.testutil.assertTrue[`OK=summary`status;"summary status OK"];
.testutil.assertTrue[11=summary`steps;"summary steps = 11"];
.testutil.assertTrue[not null summary`totalPnl;"totalPnl finite"];
.testutil.assertTrue[(summary`numRebalances)>=1;"at least the initial rebalance counted"];
.testutil.assertTrue[(summary`theoreticalGammaPnlTotal)>0f;"long-gamma theoretical pnl positive in moving market"];
.testutil.assertTrue[0f=summary`txnCostTotal;"zero txnCostRate -> zero txnCostTotal"];
.testutil.assertTrue[0f=summary`financingTotal;"zero financingRate -> zero financingTotal"];

shortCfg:@[stratCfg;`optionSide;:;`short];
shortBundle:.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;shortCfg];
.testutil.assertTrue[(shortBundle[`summary;`theoreticalGammaPnlTotal])<0f;"short-gamma theoretical pnl negative in moving market"];

-1 "PASS test_strategy_gamma_scalp_synthetic: steps=",string[summary`steps],", totalPnl=",string[summary`totalPnl],", gammaTot=",string[summary`theoreticalGammaPnlTotal];
