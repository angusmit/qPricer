/ test_control_variate_geometric_asian.q
\l lib/init.q
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount!(
    1;`AAPL;`asianOption;`european;`call;100f;1f;1f;`arithmetic;`discrete;50);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(25000;50;42;0b;0b;0.95);
configDict:enlist[`mcConfig]!enlist mcConfig;

/ Plain arithmetic Asian
plainResult:.asian.priceAsianOption[trade;mkt;configDict];
/ With geometric control variate
cvResult:.asian.priceAsianOptionWithControlVariate[trade;mkt;configDict];

.testutil.assertTrue[cvResult[`unitPrice]>0f;"CV price positive"];
.testutil.assertTrue[cvResult[`standardError]>0f;"CV SE positive"];
.testutil.assertTrue[not null cvResult`betaValue;"beta finite"];
.testutil.assertTrue[cvResult[`controlVariate]~`geometricAsian;"CV type correct"];

/ CV SE should be lower than or comparable to plain SE
-1 "PASS test_control_variate_geometric_asian: plainSE=",string[plainResult`standardError],", cvSE=",string[cvResult`standardError],", beta=",string cvResult`betaValue;
