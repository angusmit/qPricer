/ test_lookback_greeks.q
\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
    1;`AAPL;`lookbackOption;`european;`call;`fixed;100f;1f;1f;50);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(50000;50;42;0b;0b;0.95);
configDict:enlist[`mcConfig]!enlist mcConfig;

callDelta:.lookback.bumpGreek[trade;mkt;configDict;`delta];
.testutil.assertTrue[callDelta>0f;"fixed call delta positive"];

callVega:.lookback.bumpGreek[trade;mkt;configDict;`vega];
.testutil.assertTrue[callVega>0f;"vega positive"];

callRho:.lookback.bumpGreek[trade;mkt;configDict;`rho];
.testutil.assertTrue[not null callRho;"rho finite"];

/ Same seed stability
callDelta2:.lookback.bumpGreek[trade;mkt;configDict;`delta];
.testutil.assertNear[callDelta;callDelta2;0.001;"delta stable with same seed"];

-1 "PASS test_lookback_greeks: delta=",string[callDelta],", vega=",string callVega;
