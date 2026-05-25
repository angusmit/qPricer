/ test_early_exercise_boundary.q - validate early exercise boundary extraction
\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`put;100f;1f;1f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

boundaryTable:.american.extractEarlyExerciseBoundary[trade;marketData;model;config];

/ 1. Table is not empty
if[0=count boundaryTable; '"FAIL: boundary table is empty"];

/ 2. Required columns present
requiredColumns:`tradeId`underlying`optionType`timePoint`remainingTime`exerciseBoundary`hasExerciseRegion;
tableColumns:cols boundaryTable;
missingColumns:requiredColumns where not requiredColumns in tableColumns;
if[0<count missingColumns; '"FAIL: missing columns: ",", " sv string missingColumns];

/ 3. At least one row has exercise region
exerciseRows:boundaryTable where boundaryTable`hasExerciseRegion;
if[0=count exerciseRows; '"FAIL: no rows with exercise region"];

/ 4. Non-null boundaries >= 0
nonNullBoundaries:(boundaryTable`exerciseBoundary) where not null boundaryTable`exerciseBoundary;
if[any nonNullBoundaries<0f; '"FAIL: negative boundary values found"];

/ 5. Non-null boundaries <= strike
strike:100f;
if[any nonNullBoundaries>strike+1.5; '"FAIL: boundary exceeds strike + spotStep"];

/ 6. Expiry row boundary close to strike
expiryRow:boundaryTable where boundaryTable[`remainingTime]<=0.001;
if[0<count expiryRow;
    expiryBoundary:(expiryRow`exerciseBoundary) 0;
    if[not null expiryBoundary;
        if[expiryBoundary<(strike-3f); '"FAIL: expiry boundary too far from strike"]]];

/ 7. Validate error for European trade
europeanTrade:@[trade;`exerciseStyle;:;`european];
europeanErr:@[{.american.extractEarlyExerciseBoundary[x;marketData;model;config];`NO_ERROR};europeanTrade;{`ERROR}];
if[europeanErr~`NO_ERROR; '"FAIL: expected error for European trade"];

/ 8. Validate error for call trade
callTrade:@[trade;`optionType;:;`call];
callTrade:@[callTrade;`exerciseStyle;:;`american];
callErr:@[{.american.extractEarlyExerciseBoundary[x;marketData;model;config];`NO_ERROR};callTrade;{`ERROR}];
if[callErr~`NO_ERROR; '"FAIL: expected error for call trade"];

exerciseRowCount:count exerciseRows;
expiryBoundaryVal:0Nf;
if[0<count expiryRow; expiryBoundaryVal:(expiryRow`exerciseBoundary) 0];

-1 "PASS test_early_exercise_boundary: rows=",string[count boundaryTable],", exerciseRows=",string[exerciseRowCount],", expiryBoundary=",string expiryBoundaryVal;
