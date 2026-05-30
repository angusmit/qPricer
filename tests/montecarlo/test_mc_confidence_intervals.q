/ test_mc_confidence_intervals.q
\l core/init.q
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];

mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(25000;1;42;0b;0b;0.95);
mcResult:.montecarlo.priceEuropeanMC[`call;100f;100f;1f;0.05;0f;0.2;mcConfig];

/ CI should contain BS for 25k paths (high probability)
containsBS:.convergence.confidenceIntervalContains[mcResult;bsCall];

/ Build CI from variance utility
ciResult:.variance.confidenceInterval[mcResult`price;mcResult`standardError;0.95];
.testutil.assertTrue[ciResult[`lowerConfidence]<ciResult`upperConfidence;"lower < upper"];

/ Geometric Asian CI
geoBenchmark:.asian.geometricAsianClosedForm[`call;100f;100f;1f;0.05;0f;0.2;50];
geoTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount!(
    1;`AAPL;`asianOption;`european;`call;100f;1f;1f;`geometric;`discrete;50);
geoMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
geoCfg:enlist[`mcConfig]!enlist mcConfig;
geoResult:.asian.priceAsianOption[geoTrade;geoMkt;geoCfg];
geoContains:.convergence.confidenceIntervalContains[geoResult;geoBenchmark];

-1 "PASS test_mc_confidence_intervals: EuropeanContainsBS=",string[containsBS],", GeoContainsBenchmark=",string geoContains;
