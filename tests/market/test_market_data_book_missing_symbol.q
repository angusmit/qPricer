/ test_market_data_book_missing_symbol.q - missing symbol produces error rows
\l lib/init.q

tradeTable:([]
    tradeId:1 2;
    underlying:`AAPL`TSLA;
    productType:`equityOption`equityOption;
    exerciseStyle:`european`european;
    optionType:`call`call;
    strike:100 200f;
    expiry:1 1f;
    notional:1000000 1000000f;
    barrierType:`none`none;
    barrierLevel:0N 0N;
    rebate:0 0f);

spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volatilityTable:([] underlying:enlist `AAPL; volatility:enlist 0.20);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
dividendTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);

marketDataBook:.marketbook.createMarketDataBook[spotTable;volatilityTable;rateTable;dividendTable];
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

portfolioPriceTable:.portfolio.priceTradesWithMarketDataBook[tradeTable;marketDataBook;model;config];
if[not 2=count portfolioPriceTable; '"FAIL: expected 2 rows"];

priceStatus:portfolioPriceTable`status;
if[not priceStatus[0]~`OK; '"FAIL: AAPL should be OK"];
if[not priceStatus[1]~`ERROR; '"FAIL: TSLA should be ERROR"];
if[not null (portfolioPriceTable`unitPrice) 1; '"FAIL: TSLA unitPrice should be null"];

portfolioScenarioTable:.portfolio.generatePortfolioScenarioReportWithMarketDataBook[tradeTable;marketDataBook;model;config];
totalScenarioRows:count portfolioScenarioTable;
if[not totalScenarioRows=10; '"FAIL: expected 10 scenario rows (9+1), got ",string totalScenarioRows];

tslaRows:portfolioScenarioTable where portfolioScenarioTable[`tradeId]=2;
if[not 1=count tslaRows; '"FAIL: TSLA should have 1 error row"];
if[not (tslaRows`status)[0]~`ERROR; '"FAIL: TSLA scenario status should be ERROR"];

-1 "PASS test_market_data_book_missing_symbol: pricingRows=",string[count portfolioPriceTable],", scenarioRows=",string totalScenarioRows;
