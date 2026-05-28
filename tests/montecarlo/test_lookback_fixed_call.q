/ test_lookback_fixed_call.q
\l lib/init.q
trade:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
    1;`AAPL;`lookbackOption;`european;`call;`fixed;100f;1f;1f;50);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(50000;50;42;0b;0b;0.95);
configDict:enlist[`mcConfig]!enlist mcConfig;

priceResult:.lookback.priceLookbackOption[trade;mkt;configDict];
.testutil.assertTrue[priceResult[`unitPrice]>0f;"fixed call price positive"];
.testutil.assertTrue[priceResult[`standardError]>0f;"SE positive"];
.testutil.assertTrue[priceResult[`lowerConfidence]<=priceResult`unitPrice;"lower CI <= price"];

/ Fixed lookback call >= European call
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
.testutil.assertTrue[priceResult[`unitPrice]>=bsCall-0.5;"lookback call >= European call"];

-1 "PASS test_lookback_fixed_call: price=",string[priceResult`unitPrice],", SE=",string priceResult`standardError;
