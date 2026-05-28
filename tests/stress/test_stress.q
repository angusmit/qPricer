/ test_stress.q - synthetic data generation and stress tests (v0.14)
\l lib/init.q

bsModel:.model.createBlackScholesModel[];
fdmConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;200;0f;1500f;`linear;1b;1b);
configDict:`model`fdmConfig`timeStepYears`bookName`valuationDate`runLabel!(
    bsModel;fdmConfig;1%252;"stressBook";2025.01.01;"stressRun");

/ 1. Generate market data book
mktBook:.stress.generateMarketDataBook[5;2025.01.01];
.testutil.assertTrue[5=count mktBook`spotTable;"book has 5 symbols"];

/ 2. Market data columns
.testutil.assertTableColumns[mktBook`spotTable;`underlying`spot;"spot columns"];
.testutil.assertTableColumns[mktBook`volatilityTable;`underlying`volatility;"vol columns"];

/ 3. Generate supported trade table
symbolList:mktBook[`spotTable]`underlying;
supportedTrades:.stress.generateSupportedTradeTable[10;symbolList;2025.01.01];
.testutil.assertTrue[10=count supportedTrades;"10 supported trades"];
.testutil.assertTrue[all supportedTrades[`exerciseStyle]=`european;"all European"];

/ 4. Supported trades price OK
supportedResult:.portfolio.priceTradesWithMarketDataBook[supportedTrades;mktBook;bsModel;fdmConfig];
supportedOk:sum supportedResult[`status]=`OK;
.testutil.assertTrue[supportedOk=count supportedTrades;"all supported trades priced OK"];

/ 5. Mixed mode has some errors
mixedTrades:.stress.generateTradeTable[20;symbolList;2025.01.01];
mixedResult:.portfolio.priceTradesWithMarketDataBook[mixedTrades;mktBook;bsModel;fdmConfig];
.testutil.assertTrue[20=count mixedResult;"mixed result has 20 rows"];

/ 6. Missing data stress
missingResult:.stress.runMissingDataStress[10;5;2;configDict];
.testutil.assertTrue[(missingResult`errorRows)>0;"missing data produces errors"];

/ 7. Bad trade injection
enrichedTable:.stress.injectBadTrades[supportedTrades;3];
.testutil.assertTrue[13=count enrichedTable;"enriched has 13 rows"];

/ 8. Supported-only portfolio stress
stressResult:.stress.runPortfolioStress[5;3;configDict;`supportedOnly];
.testutil.assertTrue[(stressResult`okRows)>0;"supported stress has OK rows"];

-1 "PASS test_stress: supportedOk=",string[supportedOk],", missingErrors=",string missingResult`errorRows;
