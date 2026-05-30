/ test_portfolio_multi_symbol_pricing.q - multi-symbol portfolio pricing
\l core/init.q

tradeTable:([]
    tradeId:1 2 3;
    underlying:`AAPL`MSFT`NVDA;
    productType:`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`european;
    optionType:`call`put`call;
    strike:100 250 800f;
    expiry:1 1 1f;
    notional:1000000 1000000 1000000f;
    barrierType:`none`none`none;
    barrierLevel:0N 0N 0N;
    rebate:0 0 0f);

spotTable:([] underlying:`AAPL`MSFT`NVDA; spot:100 250 800f);
volatilityTable:([] underlying:`AAPL`MSFT`NVDA; volatility:0.20 0.25 0.35);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
dividendTable:([] underlying:`AAPL`MSFT`NVDA; dividendYield:0 0.01 0f);

marketDataBook:.marketbook.createMarketDataBook[spotTable;volatilityTable;rateTable;dividendTable];
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;500;0f;2400f;`linear;1b;1b);

portfolioPriceTable:.portfolio.priceTradesWithMarketDataBook[tradeTable;marketDataBook;model;config];

if[not 3=count portfolioPriceTable; '"FAIL: expected 3 rows"];

statusValues:portfolioPriceTable`status;
if[not all statusValues=`OK; '"FAIL: not all status OK"];

priceValues:portfolioPriceTable`unitPrice;
if[any priceValues<=0f; '"FAIL: non-positive unitPrice"];
if[priceValues[0]=priceValues[1]; '"FAIL: AAPL and MSFT prices should differ"];

-1 "PASS test_portfolio_multi_symbol_pricing: rows=3, allStatusOK=1";
