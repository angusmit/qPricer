/ test_merton_black_scholes_limit.q - lambda=0 should match BS
\l lib/init.q
bsLimitParams:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.2;0.0;0.0;0.0;0.05;0.0);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);

mertonCall:.merton.priceEuropeanSeries[`call;100f;100f;1f;bsLimitParams;30];
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
callError:abs mertonCall-bsCall;
.testutil.assertTrue[callError<0.01;"Merton call close to BS (lambda=0)"];

mertonPut:.merton.priceEuropeanSeries[`put;100f;100f;1f;bsLimitParams;30];
bsPut:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];
putError:abs mertonPut-bsPut;
.testutil.assertTrue[putError<0.01;"Merton put close to BS (lambda=0)"];

-1 "PASS test_merton_black_scholes_limit: callError=",string[callError],", putError=",string putError;
