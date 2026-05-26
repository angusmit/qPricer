/ test_heston_black_scholes_limit.q - volOfVol=0 should match BS
\l lib/init.q
/ BS limit: xi=0, v0=theta=0.04 => constant vol sigma=0.2
bsLimitParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.0;0.0;0.05;0.0);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(25000;100;42;0b;0b;0.95);
configDict:`mcConfig`hestonParams`modelType!(mcConfig;bsLimitParams;`heston);

hestonCallResult:.heston.priceEuropean[trade;mkt;configDict];
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
callError:abs hestonCallResult[`unitPrice]-bsCall;

.testutil.assertTrue[callError<0.5;"Heston call close to BS (xi=0)"];

/ Put
tradePut:@[trade;`optionType;:;`put];
hestonPutResult:.heston.priceEuropean[tradePut;mkt;configDict];
bsPut:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];
putError:abs hestonPutResult[`unitPrice]-bsPut;

.testutil.assertTrue[putError<0.5;"Heston put close to BS (xi=0)"];

-1 "PASS test_heston_black_scholes_limit: callError=",string[callError],", putError=",string putError;
