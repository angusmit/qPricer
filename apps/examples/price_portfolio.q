/ price_portfolio.q - portfolio pricing and batch risk example
/ Usage: q examples/price_portfolio.q

\l core/init.q
-1 "qFDM v",.qfdm.version," - Portfolio Pricing Example\n";

tradeTable:([]
    tradeId:1 2 3 4;
    underlying:`AAPL`AAPL`AAPL`AAPL;
    productType:`equityOption`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`american`european;
    optionType:`call`put`put`call;
    strike:100 100 100 100f;
    expiry:1 1 1 1f;
    notional:1000000 1000000 1000000 1000000f;
    barrierType:`none`none`none`upAndOut;
    barrierLevel:0Nf 0Nf 0Nf 130f;
    rebate:0 0 0 0f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

-1 "Portfolio prices:";
portfolioPriceTable:.portfolio.priceTrades[tradeTable;marketData;model;config];
show portfolioPriceTable;
-1 "";

-1 "Portfolio Greeks (European vanilla only):";
portfolioGreeksTable:.portfolio.calculatePortfolioGreeks[tradeTable;marketData;model;config];
show portfolioGreeksTable;
-1 "";

-1 "Portfolio scenario risk (first 18 of 36 rows):";
portfolioScenarioTable:.portfolio.generatePortfolioScenarioReport[tradeTable;marketData;model;config];
show 18 sublist portfolioScenarioTable;
-1 "... (",string[count portfolioScenarioTable]," total rows)";
-1 "";

-1 "Scenario summary:";
portfolioSummaryTable:.portfolio.summarizePortfolioRisk portfolioScenarioTable;
show portfolioSummaryTable;
