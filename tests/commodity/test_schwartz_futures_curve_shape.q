\l core/init.q
/ Futures curve shape: (a) curve is positive everywhere; (b) when x0 < theta the
/ curve is monotonically increasing in T; (c) long-T asymptote is exp(theta + sigma^2/(4 kappa)).

params:`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30);
x0:log 70f;
tenorVector:0.25 0.5 1 2 5 10 20f;

curveVector:.commodity.schwartz.futuresCurve[x0;params;tenorVector];

.testutil.assertTrue[all curveVector>0f;"curve strictly positive"];
.testutil.assertTrue[all 0f<1_deltas curveVector;"monotone increasing when x0 < theta"];

kappaVal:params`meanReversionSpeed;
sigmaVal:params`volatility;
thetaVal:params`longRunLogMean;
longTAsymptote:exp thetaVal+(sigmaVal*sigmaVal)%(4f*kappaVal);
.testutil.assertNear[last curveVector;longTAsymptote;0.01;"long-T asymptote matches"];

-1 "PASS test_schwartz_futures_curve_shape: curve ",-3!curveVector,", asymptote ",string longTAsymptote;
