\l core/init.q
/ Economic-sense check (the real correctness test, not just RMSE): a BACKWARDATED
/ curve must calibrate to a net convenience yield ABOVE the risk-free rate
/ (backwardation <=> convenience yield > carry). Curves are hand-built (monotone),
/ independent of the schwartz2 model that fits them.
tenors:0.08 0.25 0.5 0.75 1.0 1.5 2.0f;
backwardatedPrices:70 68.5 66.8 65.4 64.2 62.3 60.8f;
backCurve:([] tenor:tenors; price:backwardatedPrices);
calCfg:(enlist `riskFreeRate)!enlist 0.02;
backRes:.commodity.calibrateCurve[backCurve;`schwartz2;calCfg];
.testutil.assertTrue[(backRes`netConvenienceYield)>backRes`riskFreeRate;"backwardation -> convenience yield > rate"];
.testutil.assertTrue[(backRes`impliedFrontSlope)<0f;"backwardated fit has negative front slope"];

/ Two-sided check: a CONTANGO curve calibrates to convenience yield BELOW the rate.
contangoPrices:60 60.8 61.7 62.4 63.0 64.1 65.0f;
contangoCurve:([] tenor:tenors; price:contangoPrices);
contangoRes:.commodity.calibrateCurve[contangoCurve;`schwartz2;calCfg];
.testutil.assertTrue[(contangoRes`netConvenienceYield)<contangoRes`riskFreeRate;"contango -> convenience yield < rate"];
.testutil.assertTrue[(contangoRes`impliedFrontSlope)>0f;"contango fit has positive front slope"];

-1 "PASS test_calibrate_curve_economic: backwardation cy=",(string backRes`netConvenienceYield),", contango cy=",string contangoRes`netConvenienceYield;
