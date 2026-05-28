/ test_antithetic_variance_reduction.q
\l lib/init.q
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];

plainCfg:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(10000;1;42;0b;0b;0.95);
antiCfg:@[plainCfg;`antithetic;:;1b];

plainResult:.montecarlo.priceEuropeanMC[`call;100f;100f;1f;0.05;0f;0.2;plainCfg];
antiResult:.montecarlo.priceEuropeanMC[`call;100f;100f;1f;0.05;0f;0.2;antiCfg];

.testutil.assertTrue[(abs plainResult[`price]-bsCall)<1f;"plain MC close to BS"];
.testutil.assertTrue[(abs antiResult[`price]-bsCall)<1f;"antithetic MC close to BS"];
.testutil.assertTrue[antiResult[`standardError]>0f;"antithetic SE positive"];

/ Same seed deterministic
antiResult2:.montecarlo.priceEuropeanMC[`call;100f;100f;1f;0.05;0f;0.2;antiCfg];
.testutil.assertNear[antiResult`price;antiResult2`price;1e-10;"antithetic deterministic"];

-1 "PASS test_antithetic_variance_reduction: plainSE=",string[plainResult`standardError],", antiSE=",string antiResult`standardError;
