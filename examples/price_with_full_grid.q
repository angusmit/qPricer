/ =============================================================================
/ Example: Price with full grid output
/ =============================================================================
/ Demonstrates how to retrieve and inspect the full option value surface.
/ Run from the project root: q examples/price_with_full_grid.q
/ =============================================================================

\l lib/init.q

/ Define the trade
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    3;
    `AAPL;
    `equityOption;
    `european;
    `call;
    100f;
    1f;
    1000000f
 );

/ Define flat market data
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;
    100f;
    0.05;
    0f;
    0.2
 );

/ Create model
model:.model.createBlackScholesModel[];

/ Use a coarser grid for readable output
config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot`returnFullGrid!(
    `explicit;20;100;300f;1b
 )];

/ Price with full grid
-1 "Pricing European call with full grid...";
result:.engine.priceOptionWithGrid[trade;marketData;model;config];

/ Show price result
-1 "\nPrice result:";
show result`priceResult;

/ Inspect solver result
solverResult:result`solverResult;
-1 "\nSolver metadata:";
show solverResult`metadata;

/ Show the value grid dimensions
valueGrid:solverResult`valueGrid;
-1 "\nValue grid dimensions: ",string[count valueGrid]," spot points x ",string[count first valueGrid]," time points";

/ Show option values at t=0 for a range of spot prices
-1 "\nOption values at t=0 (valuation date):";
show ([]spot:solverResult`spotGrid;optionValue:solverResult`valueAtTimeZero);

/ Show terminal payoff (last column = expiry)
-1 "\nTerminal payoff at expiry:";
show ([]spot:solverResult`spotGrid;payoff:valueGrid[;count[first valueGrid]-1]);
