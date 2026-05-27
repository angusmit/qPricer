/ test_portfolio_merton_products.q
\l lib/init.q
mertonParams:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.2;0.5;-0.1;0.3;0.05;0.0);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);

/ European call
callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;100000f);
configDict:`mertonParams`pricingMethod!(mertonParams;`series);
callResult:.merton.priceEuropean[callTrade;mkt;configDict];
.testutil.assertTrue[callResult[`status]~`OK;"call OK"];
.testutil.assertTrue[callResult[`unitPrice]>0f;"call price positive"];
.testutil.assertNear[callResult`notionalPrice;callResult[`unitPrice]*100000f;1f;"notional correct"];

/ Put
putTrade:@[callTrade;`optionType;:;`put];
putResult:.merton.priceEuropean[putTrade;mkt;configDict];
.testutil.assertTrue[putResult[`status]~`OK;"put OK"];
.testutil.assertTrue[putResult[`unitPrice]>0f;"put price positive"];

-1 "PASS test_portfolio_merton_products: callPrice=",string[callResult`unitPrice],", putPrice=",string putResult`unitPrice;
