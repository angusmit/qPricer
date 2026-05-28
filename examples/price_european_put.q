/ =============================================================================
/ Example: Price a European put option
/ =============================================================================
/ Run from the project root: q examples/price_european_put.q
/ =============================================================================

\l lib/init.q

/ Define the trade
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;
    `AAPL;
    `equityOption;
    `european;
    `put;
    100f;
    1f;
    500000f
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

/ Configure the finite-difference solver
config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f
 )];

/ Price the option
-1 "Pricing European put option...";
priceResult:.engine.priceOption[trade;marketData;model;config];
show priceResult;

-1 "";
-1 "Unit price:     ",string priceResult`unitPrice;
-1 "Notional price: ",string priceResult`notionalPrice;
