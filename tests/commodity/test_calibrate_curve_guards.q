\l core/init.q
/ Validation guards: negative/zero price (out of lognormal domain), degenerate
/ inputs, and reserved modelTypes all raise controlled errors (no crash, no fit).
calCfg:()!();
tenors:0.08 0.25 0.5 1.0f;

/ Negative-price (April-2020 WTI negative-settle regime) -> controlled error.
negCurve:([] tenor:tenors; price:50 30 -5 20f);
negErr:@[.commodity.calibrateCurve[;`schwartz2;calCfg];negCurve;{`ERROR}];
.testutil.assertTrue[negErr~`ERROR;"negative price rejected (out of lognormal domain)"];

/ Zero price -> controlled error.
zeroCurve:([] tenor:tenors; price:50 30 0 20f);
zeroErr:@[.commodity.calibrateCurve[;`schwartz2;calCfg];zeroCurve;{`ERROR}];
.testutil.assertTrue[zeroErr~`ERROR;"zero price rejected"];

/ Empty curve -> controlled error.
emptyErr:@[.commodity.calibrateCurve[;`schwartz2;calCfg];([] tenor:`float$(); price:`float$());{`ERROR}];
.testutil.assertTrue[emptyErr~`ERROR;"empty curve rejected"];

/ Single point -> controlled error.
oneErr:@[.commodity.calibrateCurve[;`schwartz2;calCfg];([] tenor:enlist 0.5; price:enlist 60f);{`ERROR}];
.testutil.assertTrue[oneErr~`ERROR;"single-point curve rejected"];

/ Two points (below the 3-parameter linear minimum) -> controlled error.
twoErr:@[.commodity.calibrateCurve[;`schwartz2;calCfg];([] tenor:0.25 0.5; price:60 59f);{`ERROR}];
.testutil.assertTrue[twoErr~`ERROR;"two-point curve rejected (need >=3)"];

/ Reserved modelType -> controlled error.
okCurve:([] tenor:0.1 0.5 1.0f; price:60 59 58f);
mrErr:@[.commodity.calibrateCurve[okCurve;;calCfg];`mrjump;{`ERROR}];
.testutil.assertTrue[mrErr~`ERROR;"mrjump modelType reserved -> controlled error"];
unkErr:@[.commodity.calibrateCurve[okCurve;;calCfg];`bogus;{`ERROR}];
.testutil.assertTrue[unkErr~`ERROR;"unknown modelType -> controlled error"];

-1 "PASS test_calibrate_curve_guards";
