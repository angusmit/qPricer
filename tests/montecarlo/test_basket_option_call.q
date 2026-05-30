/ test_basket_option_call.q
\l core/init.q

spotTable:([] underlying:`AAPL`MSFT`NVDA; spot:100 250 800f);
volTable:([] underlying:`AAPL`MSFT`NVDA; volatility:0.2 0.25 0.35);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:`AAPL`MSFT`NVDA; dividendYield:0 0.01 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
corrTable:([] sym1:`AAPL`AAPL`MSFT; sym2:`MSFT`NVDA`NVDA; correlation:0.5 0.3 0.4);

trade:`tradeId`productType`basketSymbols`basketWeights`optionType`exerciseStyle`strike`expiry`notional!(
    1;`basketOption;`AAPL`MSFT`NVDA;0.4 0.3 0.3;`call;`european;300f;1f;100000f);

mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(
    50000;50;42;0b;0b;0.95);
configDict:enlist[`mcConfig]!enlist mcConfig;

priceResult:.basket.priceBasketOption[trade;mktBook;corrTable;configDict];

.testutil.assertTrue[priceResult[`unitPrice]>0f;"basket call price positive"];
.testutil.assertTrue[priceResult[`standardError]>0f;"SE positive"];
.testutil.assertTrue[priceResult[`lowerConfidence]<=priceResult`unitPrice;"lower CI <= price"];
.testutil.assertTrue[priceResult[`upperConfidence]>=priceResult`unitPrice;"upper CI >= price"];

-1 "PASS test_basket_option_call: price=",string[priceResult`unitPrice],", SE=",string priceResult`standardError;
