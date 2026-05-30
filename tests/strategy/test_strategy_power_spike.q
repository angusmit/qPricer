\l core/init.q
/ Build an mrjump-driven futures curve with chunky jumps; powerSpike should
/ open on a spike (jump observed) or large deviation.
pathCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.05 0.1 0.25 0.5f;
    `mrjump;
    `initialLogPrice`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility`contango`riskFreeRate!(
        log 50f;2f;log 50f;0.30;30f;0.10;0.20;0f;0.02);
    8;
    1f%252f;
    11);
bundle:.strategy.path.fromFuturesCurve pathCfg;
pathTbl:bundle`frontPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `PSC_S;`PWR;`equityOption;`european;`call;50f;0.10;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `powerSpikeCapture;
stratCfg:@[stratCfg;(`curveBundle;`stepYears);:;(bundle;1f%252f)];
runBundle:.strategy.runAndSummarize[`powerSpikeCapture;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:runBundle`result;
summary:runBundle`summary;
requiredCols:`stepIndex`stepDate`spot`frontLevel`callOptionValue`delta`spikeFlag`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`stepPnl`cumulativePnl`status`message;
.testutil.assertTableColumns[resultTbl;requiredCols;"powerSpikeCapture schema"];
.testutil.assertTrue[8=count resultTbl;"8 rows"];
.testutil.assertTrue[any resultTbl`spikeFlag;"at least one spike on mrjump path"];
-1 "PASS test_strategy_power_spike: spikeCount=",string[summary`spikeCount],", premiumCaptured=",string[summary`spikePremiumCaptured];
