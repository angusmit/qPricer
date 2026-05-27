\l lib/init.q
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(10000;50;42;0b;0b;0.95);
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];

/ Limit 1: lambda=0, xi=0, v0=theta=0.04 -> BS
bsLimitParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.0;0.0;0.0;0.0;0.0;0.05;0.0);
bsLimitResult:.bates.priceEuropean[trade;mkt;`batesParams`mcConfig!(bsLimitParams;mcConfig)];
bsError:abs bsLimitResult[`unitPrice]-bsCall;
.testutil.assertTrue[bsError<1.0;"Bates BS limit within 1.0"];

/ Limit 2: lambda=0 -> Heston-like
hestonLimitParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;0.0;0.0;0.0;0.05;0.0);
hestonLimitResult:.bates.priceEuropean[trade;mkt;`batesParams`mcConfig!(hestonLimitParams;mcConfig)];
.testutil.assertTrue[hestonLimitResult[`unitPrice]>0f;"Heston limit positive"];

-1 "PASS test_bates_limit_cases: bsError=",string[bsError],", hestonLimit=",string hestonLimitResult`unitPrice;
