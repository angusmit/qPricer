\l lib/init.q
pathCfg:`names`weights`spot0`vols`correlationMatrix`steps`stepYears`riskFreeRate`dividendYields`seed!(
    `A`B;0.5 0.5;100 100f;0.20 0.30;(1 0.3f;0.3 1f);6;1f%252f;0.02;0 0f;42);
bundle:.strategy.path.fromCorrelated pathCfg;
pathTbl:bundle`indexPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `DS_S;`IDX;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `dispersion;
stratCfg:@[stratCfg;(`pathBundle;`indexVol;`stepYears);:;(bundle;0.30;1f%252f)];
runBundle:.strategy.runAndSummarize[`dispersion;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:runBundle`result;
summary:runBundle`summary;
requiredCols:`stepIndex`stepDate`spot`indexStraddleValue`constituentStraddleSum`positionValue`netDelta`hedgeNotional`txnCost`positionPnl`hedgePnl`financingPnl`stepPnl`cumulativePnl`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"dispersion schema"];
.testutil.assertTrue[6=count resultTbl;"6 rows"];
.testutil.assertTrue[2=summary`numConstituents;"2 constituents"];
.testutil.assertTrue[all 0f<resultTbl`indexStraddleValue;"index straddle positive"];
.testutil.assertTrue[all 0f<resultTbl`constituentStraddleSum;"constituent straddle sum positive"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];
.testutil.assertTrue[(summary`direction)=`shortIndex;"direction default shortIndex"];
-1 "PASS test_strategy_dispersion: impliedCorr=",string[summary`impliedCorrAtEntry],", realizedCorr=",string[summary`realizedCorr];
