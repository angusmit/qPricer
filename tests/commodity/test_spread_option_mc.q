\l lib/init.q
/ Monte Carlo spread option vs Kirk closed form. Two correlated lognormal
/ forwards, terminal payoff max(F1_T-F2_T-K,0). F1=75,F2=70,K=3,T=0.25,
/ v1=0.30,v2=0.25,rho=0.8,r=0.05.
spreadParams:`fwd1`fwd2`strike`expiry`vol1`vol2`correlation`riskFreeRate!(75f;70f;3f;0.25;0.30;0.25;0.8;0.05);
kirkCall:.commodity.spread.kirkPrice[`call;spreadParams];
mcConfig:`pathCount`randomSeed`antithetic`confidenceLevel!(60000;42;1b;0.95);
mcResult:.commodity.spread.spreadOptionMC[`call;spreadParams;mcConfig];
mcCall:mcResult`price;
.testutil.assertTrue[mcCall>0f;"MC spread call positive"];
.testutil.assertTrue[(mcResult`pathCount)>=60000;"MC pathCount reported"];
.testutil.assertTrue[(mcResult`standardError)>0f;"MC standard error positive"];
/ Kirk is an approximation; MC is the unbiased reference. They should agree to a
/ few standard errors. Use 4*SE as the band (very safe with 60k antithetic paths).
band:4f*mcResult`standardError;
.testutil.assertTrue[(abs mcCall-kirkCall)<band;"MC call within 4 SE of Kirk"];

/ At strike 0 the MC must agree with the exact Margrabe closed form.
spreadParams0:@[spreadParams;`strike;:;0f];
margrabeParams:`fwd1`fwd2`expiry`vol1`vol2`correlation`riskFreeRate!(75f;70f;0.25;0.30;0.25;0.8;0.05);
margrabeCall:.commodity.spread.margrabePrice[`call;margrabeParams];
mcResult0:.commodity.spread.spreadOptionMC[`call;spreadParams0;mcConfig];
.testutil.assertTrue[(abs (mcResult0`price)-margrabeCall)<4f*mcResult0`standardError;"MC@0 within 4 SE of Margrabe"];

/ Determinism: same seed -> identical price.
mcRepeat:.commodity.spread.spreadOptionMC[`call;spreadParams;mcConfig];
.testutil.assertNear[mcCall;mcRepeat`price;1e-12;"MC deterministic for fixed seed"];

/ Put MC also positive; MC forward parity call-put ~ disc*(F1-F2-K).
mcPut:(.commodity.spread.spreadOptionMC[`put;spreadParams;mcConfig])`price;
.testutil.assertTrue[mcPut>0f;"MC spread put positive"];

-1 "PASS test_spread_option_mc: mcCall=",(string mcCall),", kirk=",(string kirkCall),", se=",string mcResult`standardError;
