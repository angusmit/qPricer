\l lib/init.q
/ mrjump validation rejects malformed inputs.

goodParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0f;0.25);

badParams1:@[goodParams;`meanReversionSpeed;:;-2f];
r1:@[.commodity.mrjump.validateParams;badParams1;{`ERROR}];
.testutil.assertTrue[r1~`ERROR;"negative meanReversionSpeed rejected"];

badParams2:@[goodParams;`meanReversionSpeed;:;0f];
r2:@[.commodity.mrjump.validateParams;badParams2;{`ERROR}];
.testutil.assertTrue[r2~`ERROR;"zero meanReversionSpeed rejected"];

badParams3:@[goodParams;`volatility;:;-0.1];
r3:@[.commodity.mrjump.validateParams;badParams3;{`ERROR}];
.testutil.assertTrue[r3~`ERROR;"negative volatility rejected"];

badParams4:@[goodParams;`jumpIntensity;:;-0.1];
r4:@[.commodity.mrjump.validateParams;badParams4;{`ERROR}];
.testutil.assertTrue[r4~`ERROR;"negative jumpIntensity rejected"];

badParams5:@[goodParams;`jumpVolatility;:;-0.1];
r5:@[.commodity.mrjump.validateParams;badParams5;{`ERROR}];
.testutil.assertTrue[r5~`ERROR;"negative jumpVolatility rejected"];

requiredKeys:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility;
{[missingKey;allParams]
    partialParams:(key allParams)except missingKey;
    partialDict:partialParams!allParams partialParams;
    r:@[.commodity.mrjump.validateParams;partialDict;{`ERROR}];
    .testutil.assertTrue[r~`ERROR;"missing ",string[missingKey]," rejected"];
 }[;goodParams] each requiredKeys;

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;500];
mcConfig:@[mcConfig;`timeStepCount;:;10];

r6:@[{[p] .commodity.mrjump.europeanOptionPriceMC[`straddle;log 70f;75f;1f;0.05;p;mcConfig]};goodParams;{`ERROR}];
.testutil.assertTrue[r6~`ERROR;"bad option type rejected"];

r7:@[{[p] .commodity.mrjump.europeanOptionPriceMC[`call;log 70f;-1f;1f;0.05;p;mcConfig]};goodParams;{`ERROR}];
.testutil.assertTrue[r7~`ERROR;"non-positive strike rejected"];

r8:@[{[p] .commodity.mrjump.simulatePaths[log 70f;p;-1f;mcConfig]};goodParams;{`ERROR}];
.testutil.assertTrue[r8~`ERROR;"non-positive expiry rejected"];

-1 "PASS test_mrjump_validation";
