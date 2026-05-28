\l lib/init.q
/ Schwartz validation rejects malformed inputs.

badParams1:`meanReversionSpeed`longRunLogMean`volatility!(-1.5;log 75f;0.30);
r1:@[.commodity.schwartz.validateParams;badParams1;{`ERROR}];
.testutil.assertTrue[r1~`ERROR;"negative kappa rejected"];

badParams2:`meanReversionSpeed`longRunLogMean`volatility!(0f;log 75f;0.30);
r2:@[.commodity.schwartz.validateParams;badParams2;{`ERROR}];
.testutil.assertTrue[r2~`ERROR;"zero kappa rejected"];

badParams3:`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;-0.30);
r3:@[.commodity.schwartz.validateParams;badParams3;{`ERROR}];
.testutil.assertTrue[r3~`ERROR;"negative volatility rejected"];

badParams4:`meanReversionSpeed`volatility!(1.5;0.30);
r4:@[.commodity.schwartz.validateParams;badParams4;{`ERROR}];
.testutil.assertTrue[r4~`ERROR;"missing longRunLogMean key rejected"];

goodParams:`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30);
r5:@[{[p] .commodity.schwartz.europeanOptionPrice[`straddle;log 75f;75f;1f;0.05;p]};goodParams;{`ERROR}];
.testutil.assertTrue[r5~`ERROR;"bad option type rejected"];

r6:@[{[p] .commodity.schwartz.europeanOptionPrice[`call;log 75f;-1f;1f;0.05;p]};goodParams;{`ERROR}];
.testutil.assertTrue[r6~`ERROR;"non-positive strike rejected"];

r7:@[{[p] .commodity.schwartz.europeanOptionPrice[`call;log 75f;75f;-1f;0.05;p]};goodParams;{`ERROR}];
.testutil.assertTrue[r7~`ERROR;"non-positive expiry rejected"];

-1 "PASS test_schwartz_validation";
