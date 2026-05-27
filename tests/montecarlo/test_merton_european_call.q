/ test_merton_european_call.q
\l lib/init.q
mertonParams:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.2;0.5;-0.1;0.3;0.05;0.0);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
configDict:`mertonParams`pricingMethod!(mertonParams;`series);

priceResult:.merton.priceEuropean[trade;mkt;configDict];
.testutil.assertTrue[priceResult[`status]~`OK;"status OK"];
.testutil.assertTrue[priceResult[`unitPrice]>0f;"call price positive"];
.testutil.assertTrue[not null priceResult`unitPrice;"call price finite"];

/ Higher jump intensity should change price
highLambda:@[mertonParams;`jumpIntensity;:;2.0];
highResult:.merton.priceEuropean[trade;mkt;`mertonParams`pricingMethod!(highLambda;`series)];
.testutil.assertTrue[not priceResult[`unitPrice]=highResult`unitPrice;"jump intensity changes price"];

-1 "PASS test_merton_european_call: price=",string priceResult`unitPrice;
