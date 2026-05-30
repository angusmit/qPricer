/ test_heston_european_put.q
\l core/init.q
hestonParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;0.05;0.0);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(1;`AAPL;`equityOption;`european;`put;100f;1f;1f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(25000;100;42;0b;0b;0.95);
configDict:`mcConfig`hestonParams`modelType!(mcConfig;hestonParams;`heston);

priceResult:.heston.priceEuropean[trade;mkt;configDict];
.testutil.assertTrue[priceResult[`status]~`OK;"status OK"];
.testutil.assertTrue[priceResult[`unitPrice]>0f;"put price positive"];
.testutil.assertTrue[priceResult[`standardError]>0f;"SE positive"];

-1 "PASS test_heston_european_put: price=",string[priceResult`unitPrice],", SE=",string priceResult`standardError;
