/ smoke_test_european_call.q — quick European call validation
/ Usage: q examples/smoke_test_european_call.q

\l lib/init.q
-1 "qFDM v",.qfdm.version," — European Call Smoke Test\n";

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);

model:.model.createBlackScholesModel[];

config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

-1 "Price result:";
show .engine.priceOption[trade;marketData;model;config];

-1 "\nValidation vs Black-Scholes:";
show .validation.validateEuropeanOption[trade;marketData;model;config];
