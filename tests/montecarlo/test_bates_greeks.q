\l lib/init.q
batesParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;0.5;-0.1;0.3;0.05;0.0);
callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;25;42;0b;0b;0.95);
configDict:`batesParams`mcConfig!(batesParams;mcConfig);
callDelta:.bates.bumpGreek[callTrade;mkt;configDict;`delta];
.testutil.assertTrue[callDelta>0f;"call delta positive"];
putTrade:@[callTrade;`optionType;:;`put];
putDelta:.bates.bumpGreek[putTrade;mkt;configDict;`delta];
.testutil.assertTrue[putDelta<0f;"put delta negative"];
callVega:.bates.bumpGreek[callTrade;mkt;configDict;`vega];
.testutil.assertTrue[not null callVega;"vega finite"];
jumpSens:.bates.bumpGreek[callTrade;mkt;configDict;`jumpIntensitySensitivity];
.testutil.assertTrue[not null jumpSens;"jump sensitivity finite"];
vovSens:.bates.bumpGreek[callTrade;mkt;configDict;`volOfVolSensitivity];
.testutil.assertTrue[not null vovSens;"volOfVol sensitivity finite"];
-1 "PASS test_bates_greeks: delta=",string[callDelta],", vega=",string[callVega],", jumpSens=",string jumpSens;
