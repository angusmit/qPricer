/ test_monte_carlo_greeks.q - bump-and-reprice Greeks
\l lib/init.q

callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount!(
    1;`AAPL;`asianOption;`european;`call;100f;1f;1f;`arithmetic;`discrete;50);
putTrade:@[callTrade;`optionType;:;`put];
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
configDict:`mcConfig`model`fdmConfig!(.montecarlo.defaultMcConfig[];.model.createBlackScholesModel[];.config.defaultPricingConfig[]);

callDelta:.montecarlo.bumpGreek[callTrade;mkt;configDict;`delta];
.testutil.assertTrue[callDelta>0f;"Asian call delta positive"];

putDelta:.montecarlo.bumpGreek[putTrade;mkt;configDict;`delta];
.testutil.assertTrue[putDelta<0f;"Asian put delta negative"];

callVega:.montecarlo.bumpGreek[callTrade;mkt;configDict;`vega];
.testutil.assertTrue[callVega>0f;"Asian call vega positive"];

callRho:.montecarlo.bumpGreek[callTrade;mkt;configDict;`rho];
.testutil.assertTrue[not null callRho;"rho returns finite number"];

-1 "PASS test_monte_carlo_greeks: callDelta=",string[callDelta],", putDelta=",string putDelta;
