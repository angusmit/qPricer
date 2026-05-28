/ =============================================================================
/ Example: Price a European call option
/ =============================================================================
/ Run from the project root: q examples/price_european_call.q
/ =============================================================================

\l lib/init.q

/ Define the trade
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;
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

/ Configure the finite-difference solver
config:`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;
    200;
    2000;
    300f
 );
config:.config.createFiniteDifferenceConfig[config];

/ Price the option
-1 "Pricing European call option...";
priceResult:.engine.priceOption[trade;marketData;model;config];
show priceResult;

-1 "";
-1 "Unit price:     ",string priceResult`unitPrice;
-1 "Notional price: ",string priceResult`notionalPrice;
