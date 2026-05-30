\l core/init.q
/ Synthetic round-trip: generate a curve from KNOWN schwartz2 params via the
/ model's OWN futures function, then calibrate and assert the FREE params are
/ recovered and fitRmse ~ 0. The two sides are independent (model-generated
/ curve vs grid-search optimiser output). Tolerances are documented inline.
knownParams:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.2;0.30;0.15;-0.08;0.30);
knownShortFactor0:0.20;
knownLongFactor0:4.10;
tenors:0.08 0.25 0.5 0.75 1.0 1.5f;
marketPrices:.commodity.schwartz2.futuresCurve[knownShortFactor0;knownLongFactor0;knownParams;tenors];
marketCurve:([] tenor:tenors; price:marketPrices);

/ Fix the vols at the known values (single-snapshot under-identifies vols).
calCfg:`shortVolatility`longVolatility`correlation`riskFreeRate!(0.30;0.15;0.30;0.02);
res:.commodity.calibrateCurve[marketCurve;`schwartz2;calCfg];
cp:res`calibratedParams;

.testutil.assertTrue[(res`status)=`OK;"calibration status OK"];
/ Tolerances: free shape params to 0.02, fitRmse to 0.01 (price scale ~60-72).
.testutil.assertNear[cp`shortFactor0;knownShortFactor0;0.02;"recovered shortFactor0 (X0)"];
.testutil.assertNear[cp`meanReversionSpeed;1.2;0.02;"recovered meanReversionSpeed (kappa)"];
.testutil.assertNear[cp`longDrift;-0.08;0.02;"recovered longDrift (muY)"];
.testutil.assertNear[cp`longFactor0;knownLongFactor0;0.02;"recovered longFactor0 (Y0)"];
.testutil.assertTrue[(res`fitRmse)<0.01;"fitRmse ~ 0"];

/ Fitted curve reproduces the market curve.
fitted:(res`fittedCurve)`modelPrice;
.testutil.assertTrue[1e-3>max abs marketPrices-fitted;"fitted curve == market curve"];
/ perTenorError schema + fittedCurve schema are stable and aligned.
.testutil.assertTrue[(`tenor`marketPrice`modelPrice`error)~cols res`perTenorError;"perTenorError schema"];
.testutil.assertTrue[(count tenors)=count res`perTenorError;"perTenorError one row per tenor"];

-1 "PASS test_calibrate_curve_roundtrip: rmse=",(string res`fitRmse),", iterations=",string res`iterations;
