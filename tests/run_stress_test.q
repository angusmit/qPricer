/ run_stress_test.q - full stress test and performance diagnostics
/ Usage: q run_stress_test.q

\l lib/init.q

-1 "=============================================================================";
-1 " qFDM v",.qfdm.version," Stress Test";
-1 "=============================================================================\n";

bsModel:.model.createBlackScholesModel[];
fdmConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;200;0f;1500f;`linear;1b;1b);
configDict:`model`fdmConfig`timeStepYears`bookName`valuationDate`runLabel!(
    bsModel;fdmConfig;1%252;"stressDesk";2025.01.02;"stressRun");

/ --- Supported Portfolio Scaling ---
-1 "--- Supported Portfolio Scaling ---";
supportedScaling:.perfdiag.portfolioScalingTest[10 50 100;20;configDict;`supportedOnly];
show supportedScaling;
-1 "";

/ --- Mixed Portfolio/Error Scaling ---
-1 "--- Mixed Portfolio/Error Scaling ---";
mixedScaling:.perfdiag.portfolioScalingTest[10 50 100;20;configDict;`mixed];
show mixedScaling;
-1 "";

/ --- Scenario Scaling ---
-1 "--- Scenario Scaling ---";
scenarioScaling:.perfdiag.scenarioScalingTest[10;5 10 25 50;configDict];
show scenarioScaling;
-1 "";

/ --- Missing Market Data Stress ---
-1 "--- Missing Market Data Stress ---";
missingResult:.stress.runMissingDataStress[50;20;5;configDict];
-1 "tradeCount: 50, missingSymbols: 5, okRows: ",string[missingResult`okRows],", errorRows: ",string missingResult`errorRows;
-1 "";

/ --- Grid Scaling ---
-1 "--- Grid Scaling ---";
mktBook:.stress.generateMarketDataBook[5;2025.01.01];
symbolList:mktBook[`spotTable]`underlying;
/ Use a known-good European ATM call
spotVal:(mktBook[`spotTable]`spot) 0;
gridTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;symbolList 0;`equityOption;`european;`call;spotVal;1f;1f;`none;0Nf;0f);

gridConfigs:(
    .config.mergeConfig[fdmConfig;`numberOfSpotSteps`numberOfTimeSteps!(50;100)];
    .config.mergeConfig[fdmConfig;`numberOfSpotSteps`numberOfTimeSteps!(100;500)];
    .config.mergeConfig[fdmConfig;`numberOfSpotSteps`numberOfTimeSteps!(200;1000)];
    .config.mergeConfig[fdmConfig;`numberOfSpotSteps`numberOfTimeSteps!(400;2000)]);
gridResult:.perfdiag.gridScalingTest[gridTrade;mktBook;bsModel;gridConfigs];
show gridResult;
-1 "";

/ --- Slow Step Diagnostics ---
-1 "--- Slow Step Diagnostics ---";
tradeTable:.stress.generateSupportedTradeTable[20;symbolList;2025.01.01];
prevBook:.stress.generateMarketDataBook[5;2024.12.31];
batchTimings:.perfdiag.timeBatchRun[tradeTable;mktBook;prevBook;configDict];
show batchTimings;
-1 "";
timingSummary:.perfdiag.summariseTimings batchTimings;
-1 "Total batch time: ",string[timingSummary`totalMs],"ms";

slowSteps:.perfdiag.detectSlowSteps[batchTimings;1000];
if[0<count slowSteps;
    -1 "\nSlow steps (>1000ms):";
    show slowSteps];

-1 "\n=============================================================================";
-1 " Stress test completed";
-1 "=============================================================================";
