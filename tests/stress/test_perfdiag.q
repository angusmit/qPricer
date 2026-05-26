/ test_perfdiag.q - performance diagnostics
\l lib/init.q

bsModel:.model.createBlackScholesModel[];
fdmConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;200;0f;1500f;`linear;1b;1b);
configDict:`model`fdmConfig`timeStepYears`bookName`valuationDate`runLabel!(
    bsModel;fdmConfig;1%252;"perfBook";2025.01.01;"perfRun");

/ 1. timeStep returns elapsedMs
timedStep:.perfdiag.timeStep["testStep";{[argVal] argVal*argVal};enlist 5];
.testutil.assertEqual[timedStep`stepName;"testStep";"timeStep name"];
.testutil.assertTrue[timedStep[`elapsedMs]>=0;"timeStep elapsed non-negative"];

/ 2. timeBatchRun returns timing rows
mktBook:.stress.generateMarketDataBook[3;2025.01.01];
symbolList:mktBook[`spotTable]`underlying;
tradeTable:.stress.generateSupportedTradeTable[3;symbolList;2025.01.01];
prevBook:.stress.generateMarketDataBook[3;2024.12.31];
batchTimings:.perfdiag.timeBatchRun[tradeTable;mktBook;prevBook;configDict];
.testutil.assertTrue[(count batchTimings)>=4;"batch timing has 4+ rows"];

/ 3. portfolioScalingTest with supportedOnly
scalingResult:.perfdiag.portfolioScalingTest[3 5;3;configDict;`supportedOnly];
.testutil.assertEqual[count scalingResult;2;"scaling test 2 rows"];

/ 4. scenarioScalingTest includes expected/actual/dropped
scenarioScaling:.perfdiag.scenarioScalingTest[3;3 5;configDict];
.testutil.assertEqual[count scenarioScaling;2;"scenario scaling 2 rows"];
.testutil.assertTableColumns[scenarioScaling;`expectedRows`actualRows`droppedRows;"scenario columns"];
/ expectedRows = tradeCount * 9
.testutil.assertEqual[(scenarioScaling`expectedRows)0;27;"first row expectedRows = 3*9"];
/ droppedRows = expected - actual
computedDropped:(scenarioScaling`expectedRows)-(scenarioScaling`actualRows);
.testutil.assertEqual[scenarioScaling`droppedRows;computedDropped;"droppedRows consistent"];

/ 5. gridScalingTest with known-good trade
singleTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;symbolList 0;`equityOption;`european;`call;(mktBook[`spotTable]`spot)0;1f;1f;`none;0Nf;0f);
configList:(
    .config.mergeConfig[fdmConfig;`numberOfSpotSteps`numberOfTimeSteps!(50;100)];
    .config.mergeConfig[fdmConfig;`numberOfSpotSteps`numberOfTimeSteps!(100;200)]);
gridResult:.perfdiag.gridScalingTest[singleTrade;mktBook;bsModel;configList];
.testutil.assertEqual[count gridResult;2;"grid test 2 rows"];
.testutil.assertTableColumns[gridResult;`unitPrice`benchmarkPrice`absoluteError`status`errorMessage;"grid columns"];
/ Successful rows must have non-null unitPrice and absoluteError
.testutil.assertTrue[not null (gridResult`unitPrice)0;"grid unitPrice not null"];
.testutil.assertTrue[not null (gridResult`absoluteError)0;"grid absError not null"];
.testutil.assertTrue[not null (gridResult`benchmarkPrice)0;"grid benchmark not null"];
.testutil.assertEqual[(gridResult`status)0;`OK;"grid status OK"];

/ 6. detectSlowSteps and summariseTimings
timingSummary:.perfdiag.summariseTimings batchTimings;
.testutil.assertTrue[timingSummary[`totalMs]>=0;"total ms non-negative"];

-1 "PASS test_perfdiag: batchSteps=",string[count batchTimings],", gridError=",string (gridResult`absoluteError) 0;
