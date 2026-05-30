/ test_report.q - desk reporting functions (v0.14)
\l core/init.q

tradeTable:([]
    tradeId:1 2 3;
    underlying:`AAPL`AAPL`AAPL;
    productType:`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`american;
    optionType:`call`put`put;
    strike:100 100 100f;
    expiry:1 1 1f;
    notional:1000000 1000000 1000000f;
    barrierType:`none`none`none;
    barrierLevel:0N 0N 0N;
    rebate:0 0 0f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
bsModel:.model.createBlackScholesModel[];
fdmConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

/ 1. Portfolio summary
pricingResult:.portfolio.priceTrades[tradeTable;marketData;bsModel;fdmConfig];
portfolioSummary:.report.portfolioSummary pricingResult;
if[not portfolioSummary[`totalTrades]=3; '"FAIL: totalTrades should be 3"];
if[not portfolioSummary[`okTrades]=3; '"FAIL: okTrades should be 3"];
if[not portfolioSummary[`totalNotionalPrice]>0f; '"FAIL: totalNotionalPrice should be positive"];

/ 2. Risk summary (v0.14: all 3 trades get Greeks)
greekResult:.portfolio.calculatePortfolioGreeks[tradeTable;marketData;bsModel;fdmConfig];
riskSummary:.report.riskSummary greekResult;
if[not riskSummary[`totalTrades]=3; '"FAIL: risk totalTrades should be 3"];
if[not riskSummary[`okTrades]=3; '"FAIL: risk okTrades should be 3"];

/ 3. Scenario summary
scenarioResult:.portfolio.generatePortfolioScenarioReport[tradeTable;marketData;bsModel;fdmConfig];
scenarioSummary:.report.scenarioSummary scenarioResult;
if[not 9=count scenarioSummary; '"FAIL: scenario summary should have 9 rows"];

/ 4. Error summary (no errors expected)
errorRows:.report.errorSummary pricingResult;
if[not 0=count errorRows; '"FAIL: should have no pricing errors"];

/ 5. CSV export
csvPath:"/home/claude/test_export.csv";
.report.exportCsv[pricingResult;csvPath];
csvLines:read0 hsym `$csvPath;
if[not (count csvLines)>1; '"FAIL: CSV should have header + data rows"];

-1 "PASS test_report: portfolioSummary OK, riskSummary OK, scenarioSummary=",string[count scenarioSummary]," rows, csvExport OK";
