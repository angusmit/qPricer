\l lib/init.q
batesParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;0.5;-0.1;0.3;0.05;0.0);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;25;42;0b;0b;0.95);
callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(1;`AAPL;`equityOption;`european;`call;100f;1f;100000f);
configDict:`batesParams`mcConfig!(batesParams;mcConfig);
callResult:.bates.priceEuropean[callTrade;mkt;configDict];
.testutil.assertTrue[callResult[`status]~`OK;"call OK"];
.testutil.assertTrue[callResult[`unitPrice]>0f;"call price positive"];
.testutil.assertNear[callResult`notionalPrice;callResult[`unitPrice]*100000f;1f;"notional correct"];
putTrade:@[callTrade;`optionType;:;`put];
putResult:.bates.priceEuropean[putTrade;mkt;configDict];
.testutil.assertTrue[putResult[`status]~`OK;"put OK"];
.testutil.assertTrue[putResult[`unitPrice]>0f;"put price positive"];
-1 "PASS test_portfolio_bates_products: callPrice=",string[callResult`unitPrice],", putPrice=",string putResult`unitPrice;
