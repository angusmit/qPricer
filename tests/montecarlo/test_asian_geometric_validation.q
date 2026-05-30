/ test_asian_geometric_validation.q - MC geometric vs closed-form
\l core/init.q

closedFormPrice:.asian.geometricAsianClosedForm[`call;100f;100f;1f;0.05;0f;0.2;50];

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount!(
    1;`AAPL;`asianOption;`european;`call;100f;1f;1f;`geometric;`discrete;50);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(
    50000;50;42;1b;0b;0.95);
configDict:`mcConfig`model`fdmConfig!(mcConfig;.model.createBlackScholesModel[];.config.defaultPricingConfig[]);

mcResult:.asian.priceAsianOption[trade;mkt;configDict];
mcPrice:mcResult`unitPrice;
geoError:abs mcPrice-closedFormPrice;

.testutil.assertTrue[geoError<0.3;"geometric MC close to closed-form"];
.testutil.assertTrue[closedFormPrice>0f;"closed-form positive"];

-1 "PASS test_asian_geometric_validation: MC=",string[mcPrice],", closedForm=",string[closedFormPrice],", error=",string geoError;
