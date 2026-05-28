\l lib/init.q
/ With jumpIntensity = jumpMean = jumpVolatility = 0 the mrjump simulator
/ degenerates to the Schwartz one-factor OU log-price model, for which a
/ closed-form European option price exists. Verify both MC call and MC put
/ are within 4 * standardError + 1e-3 of the Schwartz closed-form.
/ Tolerance design: 4 * SE covers ~99.99% of MC noise; the additive 1e-3
/ guards against the case where SE is very small relative to discretisation
/ error in the OU exact step.

x0:log 70f;
kappaVal:1.5;
thetaVal:log 75f;
sigmaVal:0.30;
expiryVal:1f;
strikeVal:75f;
rateVal:0.05;

mrjumpParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(kappaVal;thetaVal;sigmaVal;0f;0f;0f);
schwartzParams:`meanReversionSpeed`longRunLogMean`volatility!(kappaVal;thetaVal;sigmaVal);

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;50];

callMcResult:.commodity.mrjump.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;mrjumpParams;mcConfig];
putMcResult:.commodity.mrjump.europeanOptionPriceMC[`put;x0;strikeVal;expiryVal;rateVal;mrjumpParams;mcConfig];

callClosed:.commodity.schwartz.europeanOptionPrice[`call;x0;strikeVal;expiryVal;rateVal;schwartzParams];
putClosed:.commodity.schwartz.europeanOptionPrice[`put;x0;strikeVal;expiryVal;rateVal;schwartzParams];

callDiff:abs callMcResult[`price]-callClosed;
putDiff:abs putMcResult[`price]-putClosed;
callTol:(4f*callMcResult`standardError)+1e-3;
putTol:(4f*putMcResult`standardError)+1e-3;

.testutil.assertTrue[callDiff<callTol;"mrjump MC call within 4 SE + 1e-3 of Schwartz closed-form"];
.testutil.assertTrue[putDiff<putTol;"mrjump MC put within 4 SE + 1e-3 of Schwartz closed-form"];

-1 "PASS test_mrjump_lambda_zero_price_benchmark: callDiff=",string[callDiff],", putDiff=",string[putDiff],", callSE=",string[callMcResult`standardError],", putSE=",string putMcResult`standardError;
