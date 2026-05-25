/ run_daily_pricing.q - daily pricing, risk, and PnL explain workflow
/ Usage: q examples/run_daily_pricing.q

\l lib/init.q
-1 "qFDM v",.qfdm.version," - Daily Pricing and PnL Explain\n";

/ --- Trade book ---
tradeTable:([]
    tradeId:1 2 3;
    underlying:`AAPL`AAPL`MSFT;
    productType:`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`european;
    optionType:`call`put`call;
    strike:100 100 250f;
    expiry:1 1 1f;
    notional:1000000 500000 2000000f;
    barrierType:`none`none`none;
    barrierLevel:0N 0N 0N;
    rebate:0 0 0f);

/ --- Yesterday (t0) market data ---
spotTable0:([] underlying:`AAPL`MSFT; spot:100 250f);
volTable0:([] underlying:`AAPL`MSFT; volatility:0.20 0.25);
rateTable0:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable0:([] underlying:`AAPL`MSFT; dividendYield:0 0.01f);
book0:.marketbook.createMarketDataBook[spotTable0;volTable0;rateTable0;divTable0];

/ --- Today (t1) market data: AAPL up 2%, MSFT down 1%, vol up 1pp ---
spotTable1:([] underlying:`AAPL`MSFT; spot:102 247.5f);
volTable1:([] underlying:`AAPL`MSFT; volatility:0.21 0.26);
rateTable1:([] expiry:enlist 1f; riskFreeRate:enlist 0.0505);
divTable1:([] underlying:`AAPL`MSFT; dividendYield:0 0.01f);
book1:.marketbook.createMarketDataBook[spotTable1;volTable1;rateTable1;divTable1];

bsModel:.model.createBlackScholesModel[];
cnConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;500;0f;750f;`linear;1b;1b);

/ --- Price portfolio at t1 ---
-1 "Portfolio prices (today):";
pricingResult:.portfolio.priceTradesWithMarketDataBook[tradeTable;book1;bsModel;cnConfig];
show pricingResult;
-1 "";

-1 "Portfolio summary:";
show .report.portfolioSummary pricingResult;
-1 "";

/ --- Greeks ---
-1 "Portfolio Greeks:";
greekResult:.portfolio.calculatePortfolioGreeks[tradeTable;
    `underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;102f;0.0505;0f;0.21);
    bsModel;cnConfig];
-1 "Risk summary:";
show .report.riskSummary greekResult;
-1 "";

/ --- PnL explain ---
pnlConfig:`model`fdmConfig`timeStepYears`bookName!(bsModel;cnConfig;1%252;"equityDesk");

-1 "PnL explain:";
explainResult:.pnl.explainPortfolio[tradeTable;book0;book1;pnlConfig];
show explainResult;
-1 "";

-1 "PnL aggregate:";
show .pnl.aggregateExplain explainResult;
