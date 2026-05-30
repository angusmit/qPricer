/ test_asian_arithmetic_put.q
\l core/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount!(
    1;`AAPL;`asianOption;`european;`put;100f;1f;1f;`arithmetic;`discrete;50);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(
    50000;50;42;1b;0b;0.95);
configDict:`mcConfig`model`fdmConfig!(mcConfig;.model.createBlackScholesModel[];.config.defaultPricingConfig[]);

priceResult:.asian.priceAsianOption[trade;mkt;configDict];

.testutil.assertTrue[priceResult[`unitPrice]>0f;"Asian put price positive"];
.testutil.assertTrue[priceResult[`standardError]>0f;"standard error positive"];
.testutil.assertTrue[priceResult[`lowerConfidence]<=priceResult`unitPrice;"lower CI <= price"];
.testutil.assertTrue[priceResult[`upperConfidence]>=priceResult`unitPrice;"upper CI >= price"];

-1 "PASS test_asian_arithmetic_put: price=",string[priceResult`unitPrice],", SE=",string priceResult`standardError;
