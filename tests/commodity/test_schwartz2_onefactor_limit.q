\l lib/init.q
/ With longVolatility = 0 and longDrift = 0, the two-factor model reduces to a
/ Schwartz one-factor model in log space if initial factors are mapped:
/     shortFactor0 = log(S0) - theta
/     longFactor0  = theta
/ Then schwartz2 futures price and European option price must match the v0.33
/ schwartz one-factor closed-form to floating-point precision.

S0:70f;
thetaVal:log 75f;
kappaVal:1.5;
sigmaVal:0.30;
strikeVal:72f;
expiryVal:1f;
rateVal:0.05;

params2:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(kappaVal;sigmaVal;0f;0f;0f);
shortFactor0:(log S0)-thetaVal;
longFactor0:thetaVal;

params1:`meanReversionSpeed`longRunLogMean`volatility!(kappaVal;thetaVal;sigmaVal);
x0:log S0;

fwd2:.commodity.schwartz2.futuresPrice[shortFactor0;longFactor0;params2;expiryVal];
fwd1:.commodity.schwartz.futuresPrice[x0;params1;expiryVal];
.testutil.assertNear[fwd2;fwd1;1e-10;"futures price matches v0.33 schwartz"];

call2:.commodity.schwartz2.europeanOptionPrice[`call;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;params2];
call1:.commodity.schwartz.europeanOptionPrice[`call;x0;strikeVal;expiryVal;rateVal;params1];
.testutil.assertNear[call2;call1;1e-10;"call matches v0.33 schwartz"];

put2:.commodity.schwartz2.europeanOptionPrice[`put;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;params2];
put1:.commodity.schwartz.europeanOptionPrice[`put;x0;strikeVal;expiryVal;rateVal;params1];
.testutil.assertNear[put2;put1;1e-10;"put matches v0.33 schwartz"];

-1 "PASS test_schwartz2_onefactor_limit: fwd ",string[fwd2]," vs ",string[fwd1],", call ",string[call2]," vs ",string call1;
