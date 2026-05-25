/ run_daily_pricing.q - full daily batch: pricing, risk, PnL, audit, regression
/ Usage: q examples/run_daily_pricing.q

\l lib/init.q
-1 "qFDM v",.qfdm.version," - Daily Batch Run\n";

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

/ --- Yesterday (t0) ---
spotT0:([] underlying:`AAPL`MSFT; spot:100 250f);
volT0:([] underlying:`AAPL`MSFT; volatility:0.20 0.25);
rateT0:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divT0:([] underlying:`AAPL`MSFT; dividendYield:0 0.01f);
book0:.marketbook.createMarketDataBook[spotT0;volT0;rateT0;divT0];

/ --- Today (t1): AAPL +2%, MSFT -1%, vol +1pp ---
spotT1:([] underlying:`AAPL`MSFT; spot:102 247.5f);
volT1:([] underlying:`AAPL`MSFT; volatility:0.21 0.26);
rateT1:([] expiry:enlist 1f; riskFreeRate:enlist 0.0505);
divT1:([] underlying:`AAPL`MSFT; dividendYield:0 0.01f);
book1:.marketbook.createMarketDataBook[spotT1;volT1;rateT1;divT1];

bsModel:.model.createBlackScholesModel[];
cnConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;500;0f;750f;`linear;1b;1b);

configDict:`model`fdmConfig`timeStepYears`bookName`valuationDate`runLabel!(
    bsModel;cnConfig;1%252;"equityDesk";2025.01.02;"eodRun");

/ --- Run batch ---
-1 "Running daily batch...";
runResult:.batch.runDailyPricing[tradeTable;book1;book0;configDict];

/ --- Reports ---
-1 "\nPortfolio summary:";
show runResult`portfolioSummary;

-1 "\nRisk summary:";
show runResult`riskSummary;

-1 "\nPnL explain:";
show runResult`pnlExplainResult;

-1 "\nPnL aggregate:";
show .pnl.aggregateExplain runResult`pnlExplainResult;

-1 "\nAudit record:";
show runResult`auditRecord;

-1 "\nScenario summary:";
show runResult`scenarioSummary;

/ --- Regression ---
-1 "\nRegression check (tolerance 0.5):";
previousPricing:.portfolio.priceTradesWithMarketDataBook[tradeTable;book0;bsModel;cnConfig];
regressionResult:.regression.comparePricingRuns[runResult`pricingResult;previousPricing;0.5];
show regressionResult;

/ --- Write CSVs ---
.batch.writeDailyReports[runResult;"/home/claude";"2025.01.02"];
-1 "\nDone.";
