\l lib/init.q
/ Schwartz2 validation rejects malformed inputs.

goodParams:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.35;0.15;0f;0.3);

badParams1:@[goodParams;`meanReversionSpeed;:;-1.5];
r1:@[.commodity.schwartz2.validateParams;badParams1;{`ERROR}];
.testutil.assertTrue[r1~`ERROR;"negative meanReversionSpeed rejected"];

badParams2:@[goodParams;`meanReversionSpeed;:;0f];
r2:@[.commodity.schwartz2.validateParams;badParams2;{`ERROR}];
.testutil.assertTrue[r2~`ERROR;"zero meanReversionSpeed rejected"];

badParams3:@[goodParams;`shortVolatility;:;-0.1];
r3:@[.commodity.schwartz2.validateParams;badParams3;{`ERROR}];
.testutil.assertTrue[r3~`ERROR;"negative shortVolatility rejected"];

badParams4:@[goodParams;`longVolatility;:;-0.1];
r4:@[.commodity.schwartz2.validateParams;badParams4;{`ERROR}];
.testutil.assertTrue[r4~`ERROR;"negative longVolatility rejected"];

badParams5:@[goodParams;`correlation;:;-1f];
r5:@[.commodity.schwartz2.validateParams;badParams5;{`ERROR}];
.testutil.assertTrue[r5~`ERROR;"correlation = -1 rejected"];

badParams6:@[goodParams;`correlation;:;1f];
r6:@[.commodity.schwartz2.validateParams;badParams6;{`ERROR}];
.testutil.assertTrue[r6~`ERROR;"correlation = +1 rejected"];

badParams7:`meanReversionSpeed`shortVolatility`longVolatility`longDrift!(1.5;0.35;0.15;0f);
r7:@[.commodity.schwartz2.validateParams;badParams7;{`ERROR}];
.testutil.assertTrue[r7~`ERROR;"missing correlation key rejected"];

r8:@[{[p] .commodity.schwartz2.europeanOptionPrice[`straddle;0f;log 75f;75f;1f;0.05;p]};goodParams;{`ERROR}];
.testutil.assertTrue[r8~`ERROR;"bad option type rejected"];

r9:@[{[p] .commodity.schwartz2.europeanOptionPrice[`call;0f;log 75f;-1f;1f;0.05;p]};goodParams;{`ERROR}];
.testutil.assertTrue[r9~`ERROR;"non-positive strike rejected"];

r10:@[{[p] .commodity.schwartz2.europeanOptionPrice[`call;0f;log 75f;75f;-1f;0.05;p]};goodParams;{`ERROR}];
.testutil.assertTrue[r10~`ERROR;"non-positive expiry rejected"];

-1 "PASS test_schwartz2_validation";
