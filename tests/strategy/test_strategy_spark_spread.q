\l lib/init.q
/ sparkSpread: long a power-vs-gas spark-spread option (payoff max(P - HR*G - K,0))
/ on a correlated power/gas curve bundle, delta-hedged in both futures.
pathCfg:`names`correlationMatrix`spot0s`drifts`vols`contangos`tenors`steps`stepYears`riskFreeRate`seed!(
    `power`gas;
    (1 0.4f;0.4 1f);
    50 3.5f;
    0 0f;
    0.45 0.35f;
    0.5 0.1f;
    0.08 0.25 0.5f;
    6;
    1f%252f;
    0.03;
    77);
bundle:.strategy.path.fromCorrelatedCurves pathCfg;
pathTbl:(bundle`curves)`power;
pathTbl:pathTbl`frontPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `SPK1;`POWER;`equityOption;`european;`call;0f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;40;50;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `sparkSpread;
stratCfg:@[stratCfg;(`curveBundle;`heatRate;`vol1;`vol2;`correlation;`expiry;`stepYears);:;(bundle;8f;0.45;0.35;0.4;0.25;1f%252f)];
runBundle:.strategy.runAndSummarize[`sparkSpread;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:runBundle`result;
summary:runBundle`summary;

requiredCols:`stepIndex`stepDate`spot`leg1Front`leg2Front`leg2Mult`spreadValue`optionPrice`optionValue`delta1`delta2`hedge1`hedge2`txnCost`optionPnl`hedgePnl`financingPnl`stepPnl`cumulativePnl`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"sparkSpread schema"];
.testutil.assertTrue[6=count resultTbl;"6 rows"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"all OK"];

/ Spark spread at entry = power - heatRate*gas = 50 - 8*3.5 = 22.
.testutil.assertNear[first resultTbl`spreadValue;22f;1e-9;"spreadAtEntry = power - HR*gas"];
.testutil.assertNear[summary`spreadAtEntry;22f;1e-9;"summary spreadAtEntry"];
.testutil.assertTrue[(summary`optionPriceAtEntry)>0f;"entry option price positive"];
.testutil.assertTrue[(`sparkSpread)~summary`strategyName;"strategyName reported"];

/ Hedge: at entry the option has positive delta to power and negative to gas (a
/ spread call is long the front leg, short the gas leg). Hedges are the negatives.
.testutil.assertTrue[(first resultTbl`delta1)>0f;"positive delta to power"];
.testutil.assertTrue[(first resultTbl`delta2)<0f;"negative delta to gas"];
.testutil.assertTrue[(first resultTbl`hedge1)<0f;"short power futures to hedge"];
.testutil.assertTrue[(first resultTbl`hedge2)>0f;"long gas futures to hedge"];

/ Unhedged variant runs and holds zero hedges.
unhedgedCfg:@[stratCfg;`hedgeEnabled;:;0b];
unhedgedRes:.strategy.run[`sparkSpread;trade;pathTbl;bsModel;fdmCfg;unhedgedCfg];
.testutil.assertTrue[all 0f=unhedgedRes`hedge1;"unhedged: zero power hedge"];
.testutil.assertTrue[all 0f=unhedgedRes`hedge2;"unhedged: zero gas hedge"];

-1 "PASS test_strategy_spark_spread: spreadEntry=",(string first resultTbl`spreadValue),", optPx=",(string summary`optionPriceAtEntry),", totalPnl=",string summary`totalPnl;
