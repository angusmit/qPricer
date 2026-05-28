\l lib/init.q
/ shortVariance gate open: structural checks and short-gamma sign assertions.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.30;11;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `SV_SYNTH;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;200;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `shortVariance;
stratCfg:@[stratCfg;(`forecastVol;`entryMargin;`stepYears);:;(0.10;0.02;1f%252f)];

bundle:.strategy.runAndSummarize[`shortVariance;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;

requiredCols:`stepIndex`stepDate`spot`volatility`callPrice`putPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`premiumCollected`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"shortVariance result schema"];
.testutil.assertTrue[11=count resultTbl;"result has 11 rows (one per path step)"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];
.testutil.assertTrue[all 0f<resultTbl`callPrice;"callPrice positive"];
.testutil.assertTrue[all 0f<resultTbl`putPrice;"putPrice positive"];
.testutil.assertTrue[all 0f>resultTbl`positionValue;"positionValue negative (short straddle)"];

.testutil.assertTrue[summary`gateOpen;"gate open under these params"];
.testutil.assertTrue[(summary`premiumCollected)>0f;"premium > 0"];
.testutil.assertTrue[`OK=summary`status;"summary status OK"];
.testutil.assertTrue[11=summary`steps;"steps = 11"];
.testutil.assertTrue[not null summary`totalPnl;"totalPnl finite"];
.testutil.assertTrue[(summary`theoreticalGammaPnlTotal)<0f;"short-gamma theoretical pnl negative in moving market"];
.testutil.assertTrue[(summary`thetaPnlTotal)>0f;"short-theta theta pnl positive (collecting time decay)"];

initPremium:first resultTbl`premiumCollected;
.testutil.assertNear[summary`premiumCollected;initPremium;1e-12;"summary premium matches result premium"];

-1 "PASS test_strategy_short_variance_synthetic: premium=",string[summary`premiumCollected],", gammaTot=",string[summary`theoreticalGammaPnlTotal],", thetaTot=",string[summary`thetaPnlTotal];
