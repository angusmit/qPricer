\l core/init.q
/ Seasonality monthly-factor FIT recovers an injected monthly pattern (known-
/ answer): build a multi-year series whose value is a known function of the
/ calendar month, fit, and assert the 12 recovered factors match. Independent of
/ the forward overlay (.commodity.seasonality.factor is unchanged).
knownMonthly:0.10 0.08 0.04 0.0 -0.04 -0.08 -0.10 -0.08 -0.04 0.0 0.04 0.08;
/ 3 years of monthly observations; timeYears fractional part selects the month.
timeYears:(til 36)%12f;
monthIdx:floor 12f*timeYears-floor timeYears;
values:knownMonthly monthIdx;
fitted:.commodity.seasonality.fitMonthlyFactors[timeYears;values];
.testutil.assertTrue[12=count fitted;"12 monthly factors"];
.testutil.assertTrue[1e-12>max abs fitted-knownMonthly;"fit recovers the injected monthly pattern"];

/ The fitted factors are immediately usable by the overlay: factor at an integer
/ year (month 0) equals the January factor.
mcfg:`amplitude`phaseYears`monthlyFactors!(0f;0f;fitted);
.testutil.assertNear[.commodity.seasonality.factor[0f;mcfg];knownMonthly 0;1e-12;"fitted factors feed the overlay"];

/ Months with no observation -> null (e.g. only a half-year of data).
sparseFit:.commodity.seasonality.fitMonthlyFactors[0 0.04 0.08f;1.0 2.0 3.0];
.testutil.assertTrue[null sparseFit 6;"month with no observation is null"];

/ Length mismatch is a controlled error.
e:@[.commodity.seasonality.fitMonthlyFactors[0 0.5f;];1.0 2.0 3.0;{`ERROR}];
.testutil.assertTrue[e~`ERROR;"length mismatch rejected"];

-1 "PASS test_seasonality_fit";
