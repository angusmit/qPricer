/ test_lookback_fixed_put.q
\l core/init.q
trade:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
    1;`AAPL;`lookbackOption;`european;`put;`fixed;100f;1f;1f;50);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(50000;50;42;0b;0b;0.95);
configDict:enlist[`mcConfig]!enlist mcConfig;

priceResult:.lookback.priceLookbackOption[trade;mkt;configDict];
.testutil.assertTrue[priceResult[`unitPrice]>0f;"fixed put price positive"];
.testutil.assertTrue[priceResult[`standardError]>0f;"SE positive"];

bsPut:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];
.testutil.assertTrue[priceResult[`unitPrice]>=bsPut-0.5;"lookback put >= European put"];

-1 "PASS test_lookback_fixed_put: price=",string[priceResult`unitPrice],", SE=",string priceResult`standardError;
