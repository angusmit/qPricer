/ test_crank_nicolson_vs_explicit.q - compare CN and explicit methods
\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];

explicitCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
cnCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;500;0f;300f;`linear;1b;1b);

explicitResult:.engine.priceOption[trade;marketData;model;explicitCfg];
cnResult:.engine.priceOption[trade;marketData;model;cnCfg];
bsPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];

explicitPrice:explicitResult`unitPrice;
cnPrice:cnResult`unitPrice;
methodDiff:abs explicitPrice-cnPrice;

/ Both close to BS
if[(abs explicitPrice-bsPrice)>0.05; '"FAIL: explicit too far from BS"];
if[(abs cnPrice-bsPrice)>0.05; '"FAIL: CN too far from BS"];
/ Methods close to each other
if[methodDiff>0.10; '"FAIL: explicit vs CN diff too large: ",string methodDiff];

-1 "PASS test_crank_nicolson_vs_explicit: explicit=",string[explicitPrice],", CN=",string[cnPrice],", BS=",string[bsPrice],", diff=",string methodDiff;
