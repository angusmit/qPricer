\l core/init.q
/ Seasonality factor known-answers. Smooth cosine cycle: factor = amplitude *
/ cos(2*pi*(timeYears - phaseYears)). At phase 0 the peak is at integer years.
cfg:`amplitude`phaseYears!(0.1;0f);
.testutil.assertNear[.commodity.seasonality.factor[0f;cfg];0.1;1e-12;"t=0 -> +amplitude (winter peak)"];
.testutil.assertNear[.commodity.seasonality.factor[0.5;cfg];-0.1;1e-12;"t=0.5 -> -amplitude (summer trough)"];
.testutil.assertNear[.commodity.seasonality.factor[1f;cfg];0.1;1e-12;"t=1 -> annual cycle repeats"];

/ Phase shift moves the peak to phaseYears.
cfg2:`amplitude`phaseYears!(0.2;0.25);
.testutil.assertNear[.commodity.seasonality.factor[0.25;cfg2];0.2;1e-12;"peak shifts to phaseYears"];

/ Vectorised over a tenor vector.
fv:.commodity.seasonality.factor[0 0.5 1f;cfg];
.testutil.assertTrue[3=count fv;"vectorised over tenor vector"];
.testutil.assertNear[fv 1;-0.1;1e-12;"vector element t=0.5"];

/ Explicit 12-month override takes precedence over the cosine cycle.
monthly:0.05 0.04 0.03 0.0 -0.02 -0.04 -0.05 -0.04 -0.02 0.0 0.02 0.04;
mcfg:`amplitude`phaseYears`monthlyFactors!(0f;0f;monthly);
.testutil.assertNear[.commodity.seasonality.factor[0f;mcfg];0.05;1e-12;"month 0 (Jan) factor"];
.testutil.assertNear[.commodity.seasonality.factor[0.5f;mcfg];-0.05;1e-12;"month 6 (Jul) factor"];
.testutil.assertNear[.commodity.seasonality.factor[1.0f;mcfg];0.05;1e-12;"wraps to month 0 at t=1"];

/ A non-12 monthlyFactors vector is rejected.
badCfg:`amplitude`phaseYears`monthlyFactors!(0f;0f;1 2 3f);
bad:@[.commodity.seasonality.factor[0f;];badCfg;{`ERROR}];
.testutil.assertTrue[bad~`ERROR;"non-12 monthlyFactors rejected"];

-1 "PASS test_seasonality_factor";
