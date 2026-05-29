\l lib/init.q
/ crackSpread: long a refined-product-vs-crude crack-spread option (payoff
/ max(P_product - crackRatio*P_crude - K, 0)) on a correlated product/crude
/ curve bundle, delta-hedged in both futures. 1:1 crack here (crackRatio=1).
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
    `CRK1;`PRODUCT;`equityOption;`european;`call;0f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `crackSpread;
stratCfg:@[stratCfg;(`curveBundle;`crackRatio;`vol1;`vol2;`correlation;`expiry;`stepYears);:;(bundle;1f;0.35;0.30;0.6;0.25;1f%252f)];
runBundle:.strategy.runAndSummarize[`crackSpread;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:runBundle`result;
summary:runBundle`summary;

requiredCols:`stepIndex`stepDate`spot`leg1Front`leg2Front`leg2Mult`spreadValue`optionPrice`optionValue`delta1`delta2`hedge1`hedge2`txnCost`optionPnl`hedgePnl`financingPnl`stepPnl`cumulativePnl`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"crackSpread schema"];
.testutil.assertTrue[6=count resultTbl;"6 rows"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];

/ Crack spread at entry = product - crackRatio*crude = 95 - 1*80 = 15.
.testutil.assertNear[first resultTbl`spreadValue;15f;1e-9;"spreadAtEntry = product - ratio*crude"];
.testutil.assertNear[summary`spreadAtEntry;15f;1e-9;"summary spreadAtEntry"];
.testutil.assertTrue[(summary`optionPriceAtEntry)>0f;"entry option price positive"];
.testutil.assertTrue[(`crackSpread)~summary`strategyName;"strategyName reported"];
.testutil.assertNear[first resultTbl`leg2Mult;1f;1e-12;"crackRatio reported as leg2Mult"];

/ A larger crack ratio (e.g. 1.2) lowers the spread and thus the call value.
wideCfg:@[stratCfg;`crackRatio;:;1.2];
wideRes:.strategy.runAndSummarize[`crackSpread;trade;pathTbl;bsModel;fdmCfg;wideCfg];
.testutil.assertTrue[(wideRes[`summary;`spreadAtEntry])<summary`spreadAtEntry;"higher crackRatio -> lower spread"];

-1 "PASS test_strategy_crack_spread: spreadEntry=",(string first resultTbl`spreadValue),", optPx=",(string summary`optionPriceAtEntry),", totalPnl=",string summary`totalPnl;
