/ test_benchmark.q - benchmark module sanity tests
\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];

/ --- estimateGridMemory ---
testCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;100;1000;0f;300f;`linear;1b;1b);
memEst:.benchmark.estimateGridMemory testCfg;
expectedGridPoints:101*1001;
if[not memEst[`gridPointCount]=expectedGridPoints;
    '"FAIL: gridPointCount expected ",string[expectedGridPoints]," got ",string memEst`gridPointCount];
if[not memEst[`estimatedGridMegabytes]>0f;
    '"FAIL: estimatedGridMegabytes not positive"];

/ --- runPricingBenchmark ---
configList:(
    `method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
        `explicit;100;1000;0f;300f;`linear;1b;1b);
    `method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
        `crankNicolson;100;250;0f;300f;`linear;1b;1b));

benchmarkTable:.benchmark.runPricingBenchmark[trade;marketData;model;configList];
rowCount:count benchmarkTable;

/ 1. Correct row count
if[not rowCount=2; '"FAIL: expected 2 rows, got ",string rowCount];

/ 2. Required columns
requiredCols:`method`numberOfSpotSteps`numberOfTimeSteps`gridPointCount`estimatedGridMegabytes`runtimeMilliseconds`unitPrice`closedFormPrice`absoluteError`relativeError;
tableCols:cols benchmarkTable;
missingCols:requiredCols where not requiredCols in tableCols;
if[0<count missingCols; '"FAIL: missing columns: ",", " sv string missingCols];

/ 3. Runtime non-negative
runtimeValues:benchmarkTable`runtimeMilliseconds;
if[any runtimeValues<0; '"FAIL: negative runtime"];

/ 4. Prices positive
if[any (benchmarkTable`unitPrice)<=0f; '"FAIL: non-positive unitPrice"];
if[any (benchmarkTable`closedFormPrice)<=0f; '"FAIL: non-positive closedFormPrice"];

/ 5. Absolute error reasonable
if[any (benchmarkTable`absoluteError)>0.25; '"FAIL: absolute error too large"];

/ 6. Methods match input
expectedMethods:`explicit`crankNicolson;
if[not (benchmarkTable`method)~expectedMethods; '"FAIL: methods mismatch"];

explicitError:(benchmarkTable`absoluteError) 0;
cnError:(benchmarkTable`absoluteError) 1;

-1 "PASS test_benchmark: rows=",string[rowCount],", explicitError=",string[explicitError],", crankNicolsonError=",string cnError;
