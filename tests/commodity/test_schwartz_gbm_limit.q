\l core/init.q
/ As kappa -> 0 with longRunLogMean = log S0, Schwartz reduces to a lognormal
/ terminal price with mean log S0 and variance sigma^2 T. The European option
/ price then equals Black-76 on forward S0 * exp(sigma^2 T / 2) at the same vol.

S0:75f;
strikeVal:75f;
expiryVal:1f;
sigmaVal:0.30;
rateVal:0f;
kappaSmall:1e-7;
x0:log S0;
params:`meanReversionSpeed`longRunLogMean`volatility!(kappaSmall;x0;sigmaVal);

schwartzCall:.commodity.schwartz.europeanOptionPrice[`call;x0;strikeVal;expiryVal;rateVal;params];
schwartzPut:.commodity.schwartz.europeanOptionPrice[`put;x0;strikeVal;expiryVal;rateVal;params];

implicitForward:S0*exp 0.5*sigmaVal*sigmaVal*expiryVal;
b76Call:.commodity.black76.price[`call;implicitForward;strikeVal;expiryVal;sigmaVal;rateVal];
b76Put:.commodity.black76.price[`put;implicitForward;strikeVal;expiryVal;sigmaVal;rateVal];

.testutil.assertNear[schwartzCall;b76Call;1e-3;"GBM limit call matches Black-76"];
.testutil.assertNear[schwartzPut;b76Put;1e-3;"GBM limit put matches Black-76"];

-1 "PASS test_schwartz_gbm_limit: call ",string[schwartzCall]," vs ",string[b76Call],", put ",string[schwartzPut]," vs ",string b76Put;
