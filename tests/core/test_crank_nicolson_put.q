/ test_crank_nicolson_put.q - validate CN put price vs Black-Scholes
\l core/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
cnConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;500;0f;300f;`linear;1b;1b);

cnResult:.engine.priceOption[trade;marketData;model;cnConfig];
cnUnitPrice:cnResult`unitPrice;
bsPrice:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];
absError:abs cnUnitPrice-bsPrice;
tolerance:0.05;

if[absError>tolerance; '"FAIL test_crank_nicolson_put: error=",string[absError]," > tolerance=",string tolerance];
-1 "PASS test_crank_nicolson_put: CN=",string[cnUnitPrice],", BS=",string[bsPrice],", error=",string absError;
