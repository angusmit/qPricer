\l core/init.q
/ jumpPremium with Merton (series price, fast) at entry, gate open since
/ a chunky jump intensity in the jumpModel will push its price away from BS.
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.20;5;1f%252f;0.02;0f;42);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `JP_S;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `jumpPremium;
stratCfg:@[stratCfg;`stepYears;:;1f%252f];
/ Crank up jumpIntensity so jumpPx > bsPx significantly, threshold low enough that gate opens
stratCfg:@[stratCfg;`jumpModelParams;:;
    `volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(
        0.20;3f;-0.05;0.30;0.02;0f)];
stratCfg:@[stratCfg;`premiumThreshold;:;0.10];
bundle:.strategy.runAndSummarize[`jumpPremium;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;
requiredCols:`stepIndex`stepDate`spot`volatility`optionPrice`delta`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`bsPrice0`jumpModelPrice0`jumpPremium`gateOpen`tradeSide`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"jumpPremium schema"];
.testutil.assertTrue[5=count resultTbl;"5 rows"];
.testutil.assertTrue[summary`gateOpen;"gate open: |premium| > threshold"];
.testutil.assertTrue[(summary`jumpModelPrice0)>summary`bsPrice0;"jumpPx0 > bsPx0 (jump intensity raises premium)"];
.testutil.assertTrue[(summary`tradeSide)=`short;"auto direction sells when jump model prices above BS"];
.testutil.assertTrue[(summary`jumpModelName)=`merton;"jumpModelName=merton"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];
-1 "PASS test_strategy_jump_premium: bsPx0=",string[summary`bsPrice0],", jumpPx0=",string[summary`jumpModelPrice0],", premium=",string[summary`jumpPremium];
