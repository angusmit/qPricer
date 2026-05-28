\l lib/init.q
/ For an OU log-price with no jumps, transitionMomentsApprox should:
/   (a) at large tau, return variance ~ sigma^2 / (2 * kappa) and mean ~ theta
/   (b) at moderate tau, return variance = sigma^2 * (1 - exp(-2 kappa tau)) / (2 kappa)
/ Both branches reduce mrjump moments to the closed-form OU integrated variance.

kappaVal:1.5;
thetaVal:log 75f;
sigmaVal:0.30;
x0:log 50f;
largeTauVal:50f;
moderateTauVal:1f;

zeroJumpParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(kappaVal;thetaVal;sigmaVal;0f;0f;0f);

targetStationaryVariance:(sigmaVal*sigmaVal)%(2f*kappaVal);

largeMoments:.commodity.mrjump.transitionMomentsApprox[x0;zeroJumpParams;largeTauVal];
.testutil.assertNear[largeMoments`variance;targetStationaryVariance;1e-10;"long-time variance approaches sigma^2/(2 kappa)"];
.testutil.assertNear[largeMoments`mean;thetaVal;1e-10;"long-time mean approaches theta"];

decay2Val:exp neg 2f*kappaVal*moderateTauVal;
expectedModerateVariance:(sigmaVal*sigmaVal)*(1f-decay2Val)%(2f*kappaVal);
moderateMoments:.commodity.mrjump.transitionMomentsApprox[x0;zeroJumpParams;moderateTauVal];
.testutil.assertNear[moderateMoments`variance;expectedModerateVariance;1e-12;"moderate-tau variance matches closed-form OU integrated variance"];

decayMean:exp neg kappaVal*moderateTauVal;
expectedModerateMean:(decayMean*x0)+(1f-decayMean)*thetaVal;
.testutil.assertNear[moderateMoments`mean;expectedModerateMean;1e-12;"moderate-tau mean matches closed-form OU conditional mean"];

-1 "PASS test_mrjump_stationary_variance: largeTauVariance=",string[largeMoments`variance],", target=",string[targetStationaryVariance],", largeTauMean=",string[largeMoments`mean],", theta=",string thetaVal;
