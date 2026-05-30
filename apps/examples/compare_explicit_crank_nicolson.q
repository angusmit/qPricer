/ compare_explicit_crank_nicolson.q - compare explicit and CN solvers
/ Usage: q examples/compare_explicit_crank_nicolson.q

\l core/init.q
-1 "qFDM v",.qfdm.version," - Explicit vs Crank-Nicolson\n";

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];

explicitCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
cnCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;500;0f;300f;`linear;1b;1b);

bsPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
explicitPrice:(.engine.priceOption[trade;marketData;model;explicitCfg])`unitPrice;
cnPrice:(.engine.priceOption[trade;marketData;model;cnCfg])`unitPrice;

-1 "European call (S=100, K=100, T=1, r=5%, vol=20%):";
-1 "  Black-Scholes:   ",string bsPrice;
-1 "  Explicit (200x2000): ",string[explicitPrice],"  error=",string abs explicitPrice-bsPrice;
-1 "  CN (200x500):        ",string[cnPrice],"  error=",string abs cnPrice-bsPrice;
-1 "";

/ Now put
putTrade:@[trade;`optionType;:;`put];
putTrade:@[putTrade;`tradeId;:;2];
bsPut:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];
explicitPut:(.engine.priceOption[putTrade;marketData;model;explicitCfg])`unitPrice;
cnPut:(.engine.priceOption[putTrade;marketData;model;cnCfg])`unitPrice;

-1 "European put:";
-1 "  Black-Scholes:   ",string bsPut;
-1 "  Explicit (200x2000): ",string[explicitPut],"  error=",string abs explicitPut-bsPut;
-1 "  CN (200x500):        ",string[cnPut],"  error=",string abs cnPut-bsPut;
