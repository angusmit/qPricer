/ test_lookback_floating_call.q
\l lib/init.q
trade:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
    1;`AAPL;`lookbackOption;`european;`call;`floating;0Nf;1f;1f;50);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(50000;50;42;0b;0b;0.95);
configDict:enlist[`mcConfig]!enlist mcConfig;

priceResult:.lookback.priceLookbackOption[trade;mkt;configDict];
.testutil.assertTrue[priceResult[`unitPrice]>0f;"floating call price positive"];
.testutil.assertTrue[priceResult[`standardError]>0f;"SE positive"];
.testutil.assertTrue[priceResult[`lowerConfidence]<=priceResult`unitPrice;"lower CI <= price"];

-1 "PASS test_lookback_floating_call: price=",string[priceResult`unitPrice],", SE=",string priceResult`standardError;
