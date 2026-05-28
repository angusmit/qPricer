/ calculate_greeks.q — compute and validate Greeks
/ Usage: q examples/calculate_greeks.q

\l ./lib/init.q
-1 "qFDM v",.qfdm.version," - Greeks Example\n";

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

-1 "Price:";
show .engine.priceOption[trade;marketData;model;config];

-1 "\nFDM Greeks:";
show .greeks.calculateGreeks[trade;marketData;model;config];

-1 "\nGreeks validation (FDM vs Black-Scholes):";
show .validation.validateGreeks[trade;marketData;model;config];
