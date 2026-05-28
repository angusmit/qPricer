/ test_monte_carlo_european.q - MC European vs Black-Scholes
\l lib/init.q

mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(
    50000;1;42;0b;0b;0.95);

bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
mcCallResult:.montecarlo.priceEuropeanMC[`call;100f;100f;1f;0.05;0f;0.2;mcConfig];
mcCallPrice:mcCallResult`price;
callSE:mcCallResult`standardError;
callError:abs mcCallPrice-bsCall;

.testutil.assertTrue[callError<0.5;"MC call within 0.5 of BS"];
.testutil.assertTrue[mcCallResult[`standardError]>0f;"call SE positive"];

bsPut:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];
mcPutResult:.montecarlo.priceEuropeanMC[`put;100f;100f;1f;0.05;0f;0.2;mcConfig];
mcPutPrice:mcPutResult`price;
putError:abs mcPutPrice-bsPut;

.testutil.assertTrue[putError<0.5;"MC put within 0.5 of BS"];

-1 "PASS test_monte_carlo_european: BS call=",string[bsCall],", MC call=",string[mcCallPrice],", callError=",string[callError],", putError=",string putError;
