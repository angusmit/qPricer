\l core/init.q
/ ============================================================================
/ calibrate_crude_curve.q - real-data example (NOT a test).
/ ----------------------------------------------------------------------------
/ Loads the REAL Barchart WTI forward curve via .parser.crude, takes a
/ NORMAL-REGIME snapshot (2020-01-06, the pre-COVID backwardation curve), fits
/ the two-factor Schwartz model with .commodity.calibrateCurve, and prints the
/ recovered params, fit RMSE, fitted-vs-market curve, and the economic read.
/ ----------------------------------------------------------------------------
/ CAVEAT: the Schwartz models are log-price. The April-2020 WTI negative-settle
/ dates are OUT OF DOMAIN and would (correctly) raise a controlled error; this
/ example deliberately uses a normal-regime backwardation snapshot.
/ ----------------------------------------------------------------------------
/ The CSVs are user-supplied and gitignored - not committed. The test suite
/ exercises calibration on synthetic curves only.
/ ============================================================================
crudeDir:"data/barchart/CRUDE";
asofSnap:2020.01.06;
-1 "Loading WTI curve from ",crudeDir," , as of ",(string asofSnap)," ...";
longTable:.parser.crude.loadAll crudeDir;
curve:.parser.crude.curveAt[longTable;asofSnap];
-1 "Market curve (",(string count curve)," contracts):";
show curve;

calCfg:`shortVolatility`longVolatility`correlation`riskFreeRate!(0.35;0.20;0.30;0.02);
res:.commodity.calibrateCurve[select tenor,price from curve;`schwartz2;calCfg];

-1 "";
-1 "Calibrated schwartz2 parameters (free: X0, Y0, kappa, muY; vols fixed):";
show res`calibratedParams;
-1 "fitRmse  = ",string res`fitRmse;
-1 "iterations = ",string res`iterations;
-1 "";
-1 "Fitted vs market:";
show res`perTenorError;

cyVal:res`netConvenienceYield;
rateVal:res`riskFreeRate;
shapeRead:$[cyVal>rateVal;"EXCEEDS the rate -> BACKWARDATION (convenience yield > carry), as observed";"is below the rate -> contango"];
-1 "";
-1 "riskFreeRate        = ",string rateVal;
-1 "impliedFrontSlope   = ",string res`impliedFrontSlope;
-1 "netConvenienceYield = ",string cyVal;
-1 "Economic read: convenience yield ",shapeRead;
