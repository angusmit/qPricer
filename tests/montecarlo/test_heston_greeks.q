/ test_heston_greeks.q
\l lib/init.q
hestonParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;0.05;0.0);
callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(10000;50;42;0b;0b;0.95);
configDict:`mcConfig`hestonParams`modelType!(mcConfig;hestonParams;`heston);

callDelta:.heston.bumpGreek[callTrade;mkt;configDict;`delta];
.testutil.assertTrue[callDelta>0f;"call delta positive"];

putTrade:@[callTrade;`optionType;:;`put];
putDelta:.heston.bumpGreek[putTrade;mkt;configDict;`delta];
.testutil.assertTrue[putDelta<0f;"put delta negative"];

callVega:.heston.bumpGreek[callTrade;mkt;configDict;`vega];
.testutil.assertTrue[not null callVega;"vega finite"];

callRho:.heston.bumpGreek[callTrade;mkt;configDict;`rho];
.testutil.assertTrue[not null callRho;"rho finite"];

/ Stability
callDelta2:.heston.bumpGreek[callTrade;mkt;configDict;`delta];
.testutil.assertNear[callDelta;callDelta2;0.001;"delta stable"];

-1 "PASS test_heston_greeks: callDelta=",string[callDelta],", putDelta=",string[putDelta],", vega=",string callVega;
