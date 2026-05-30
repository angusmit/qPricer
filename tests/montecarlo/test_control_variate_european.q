/ test_control_variate_european.q - European MC with BS control variate
\l core/init.q
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];

mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(10000;1;42;0b;0b;0.95);
pathMatrix:.montecarlo.simulateGBMPaths[100f;0.05;0f;0.2;1f;mcConfig];
terminalSpots:last each pathMatrix;
callPayoff:0f|terminalSpots-100f;

/ Use discounted BS as control expected value (undiscounted)
discountFactor:exp neg 0.05;
controlExpected:bsCall%discountFactor;

betaValue:.variance.controlVariateBeta[callPayoff;callPayoff];
.testutil.assertTrue[not null betaValue;"beta is finite"];

/ Adjust: since target=control for European, CV should collapse to known value
adjustedPayoff:.variance.controlVariateAdjust[callPayoff;callPayoff;controlExpected];
adjustedPrice:discountFactor*avg adjustedPayoff;
adjustedError:abs adjustedPrice-bsCall;

.testutil.assertTrue[adjustedError<0.01;"CV European very close to BS"];

-1 "PASS test_control_variate_european: adjustedError=",string adjustedError;
