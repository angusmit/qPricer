/ price_multi_symbol_portfolio.q - multi-symbol portfolio pricing
/ Usage: q examples/price_multi_symbol_portfolio.q

\l lib/init.q
-1 "qFDM v",.qfdm.version," - Multi-Symbol Portfolio Example\n";

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
    barrierLevel:3#0Nf;
    rebate:0 0 0f);

spotTable:([] underlying:`AAPL`MSFT`NVDA; spot:100 250 800f);
volatilityTable:([] underlying:`AAPL`MSFT`NVDA; volatility:0.20 0.25 0.35);
rateTable:([] expiry:0.25 0.5 1 2f; riskFreeRate:0.045 0.0475 0.05 0.0525);
dividendTable:([] underlying:`AAPL`MSFT`NVDA; dividendYield:0 0.01 0f);

marketDataBook:.marketbook.createMarketDataBook[spotTable;volatilityTable;rateTable;dividendTable];
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;2400f;`linear;1b;0b);

-1 "Market data book:";
-1 "  Spots: AAPL=100, MSFT=250, NVDA=800";
-1 "  Vols:  AAPL=20%, MSFT=25%, NVDA=35%";
-1 "  Divs:  AAPL=0%, MSFT=1%, NVDA=0%";
-1 "";

-1 "Portfolio prices:";
portfolioPriceTable:.portfolio.priceTradesWithMarketDataBook[tradeTable;marketDataBook;model;config];
show portfolioPriceTable;
-1 "";

-1 "Portfolio scenario risk (first 9 rows):";
portfolioScenarioTable:.portfolio.generatePortfolioScenarioReportWithMarketDataBook[tradeTable;marketDataBook;model;config];
show 9 sublist portfolioScenarioTable;
-1 "... (",string[count portfolioScenarioTable]," total rows)";
