\l lib/init.q
/ Closed-gate scenario: simple curve with NO jumps and small drift -> no spike fires.
closedCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.05 0.10 0.25f;`simple;
    `spot0`drift`volatility`contango`riskFreeRate!(50f;0f;0.05;0f;0.02);
    8;1f%252f;3);
closedBundle:.strategy.path.fromFuturesCurve closedCfg;
closedPath:closedBundle`frontPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `PSC_G;`PWR;`equityOption;`european;`call;50f;0.10;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
closedStrat:.strategy.defaultConfig `powerSpikeCapture;
closedStrat:@[closedStrat;(`curveBundle;`stepYears;`deviationThreshold);:;(closedBundle;1f%252f;100f)];
closedRun:.strategy.runAndSummarize[`powerSpikeCapture;trade;closedPath;bsModel;fdmCfg;closedStrat];
.testutil.assertTrue[not closedRun[`summary;`gateOpen];"low-vol no-jump path: gate stays closed"];
.testutil.assertTrue[0=closedRun[`summary;`spikeCount];"no spikes"];
/ Open-gate: high jump intensity guarantees jumps in the curve, gate fires
openCfg:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed!(
    0.05 0.10 0.25f;`mrjump;
    `initialLogPrice`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility`contango`riskFreeRate!(
        log 50f;2f;log 50f;0.30;30f;0.10;0.20;0f;0.02);
    8;1f%252f;11);
openBundle:.strategy.path.fromFuturesCurve openCfg;
openPath:openBundle`frontPath;
openStrat:@[closedStrat;(`curveBundle;`deviationThreshold);:;(openBundle;2f)];
openRun:.strategy.runAndSummarize[`powerSpikeCapture;trade;openPath;bsModel;fdmCfg;openStrat];
.testutil.assertTrue[openRun[`summary;`gateOpen];"mrjump path: gate fires"];
.testutil.assertTrue[(openRun[`summary;`spikeCount])>0;"spike count > 0"];
-1 "PASS test_strategy_power_spike_gate: closedSpikes=0, openSpikes=",string[openRun[`summary;`spikeCount]];
