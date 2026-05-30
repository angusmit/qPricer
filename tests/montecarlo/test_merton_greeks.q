/ test_merton_greeks.q
\l core/init.q
mertonParams:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.2;0.5;-0.1;0.3;0.05;0.0);
callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
configDict:`mertonParams`pricingMethod!(mertonParams;`series);

callDelta:.merton.bumpGreek[callTrade;mkt;configDict;`delta];
.testutil.assertTrue[callDelta>0f;"call delta positive"];

putTrade:@[callTrade;`optionType;:;`put];
putDelta:.merton.bumpGreek[putTrade;mkt;configDict;`delta];
.testutil.assertTrue[putDelta<0f;"put delta negative"];

callVega:.merton.bumpGreek[callTrade;mkt;configDict;`vega];
.testutil.assertTrue[not null callVega;"vega finite"];

callRho:.merton.bumpGreek[callTrade;mkt;configDict;`rho];
.testutil.assertTrue[not null callRho;"rho finite"];

jumpSens:.merton.bumpGreek[callTrade;mkt;configDict;`jumpIntensitySensitivity];
.testutil.assertTrue[not null jumpSens;"jump intensity sensitivity finite"];

-1 "PASS test_merton_greeks: delta=",string[callDelta],", vega=",string[callVega],", jumpSens=",string jumpSens;
