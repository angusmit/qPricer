/ test_asian_arithmetic_call.q
\l core/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount!(
    1;`AAPL;`asianOption;`european;`call;100f;1f;1f;`arithmetic;`discrete;50);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(
    50000;50;42;1b;0b;0.95);
configDict:`mcConfig`model`fdmConfig!(mcConfig;.model.createBlackScholesModel[];.config.defaultPricingConfig[]);

priceResult:.asian.priceAsianOption[trade;mkt;configDict];

.testutil.assertTrue[priceResult[`unitPrice]>0f;"Asian call price positive"];
.testutil.assertTrue[priceResult[`standardError]>0f;"standard error positive"];
.testutil.assertTrue[priceResult[`lowerConfidence]<=priceResult`unitPrice;"lower CI <= price"];
.testutil.assertTrue[priceResult[`upperConfidence]>=priceResult`unitPrice;"upper CI >= price"];

bsEuropean:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
.testutil.assertTrue[priceResult[`unitPrice]<=bsEuropean+0.5;"Asian call <= European call (with tolerance)"];

-1 "PASS test_asian_arithmetic_call: price=",string[priceResult`unitPrice],", SE=",string priceResult`standardError;
