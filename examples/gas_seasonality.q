\l lib/init.q
/ ============================================================================
/ gas_seasonality.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ ----------------------------------------------------------------------------
/ IF Day_GAS_YYYYMM00.csv files are present under data/barchart/GAS, fits the
/ seasonality overlay's monthly factors to the real natural-gas forward-curve
/ seasonal pattern (winter premium / summer trough) and calibrates schwartz2 to
/ an NG snapshot. IF ABSENT, skips with a note (the overlay + fitMonthlyFactors
/ are covered by synthetic tests). CSVs are user-supplied and gitignored.
/ ============================================================================
gasDir:"data/barchart/GAS";
gasFiles:@[{key hsym `$x};gasDir;{[e] ()}];
if[0=count gasFiles;
    -1 "NG seasonality real demo SKIPPED - no data at ",gasDir," (data-conditional).";
    -1 "The seasonality overlay and .commodity.seasonality.fitMonthlyFactors are";
    -1 "exercised by tests/commodity/test_seasonality_fit.q (synthetic, known-answer).";
    -1 "Drop Day_GAS_YYYYMM00.csv files into ",gasDir," to run the real NG fit.";
    exit 0];

/ --- real branch (runs only when GAS data is present) ---
longTable:.parser.futures.loadAll[gasDir;`GAS];
allDates:asc distinct longTable`date;
curveHistory:.parser.futures.curveHistory[longTable;allDates];
/ Seasonal premium by DELIVERY month: cross-sectionally demean each date's curve
/ (log price), then average by the contract's exact delivery month (contractYM mod
/ 100). fitMonthlyFactors takes timeYears whose fractional part selects the month.
ch:select from curveHistory where price>0f;
ch:update logP:log price from ch;
ch:update demeaned:logP-(avg logP) by asofDate from ch;
deliveryMonth:(ch`contractYM) mod 100;
deliveryFracYear:(`float$deliveryMonth-1)%12f;
fitted:.commodity.seasonality.fitMonthlyFactors[deliveryFracYear;ch`demeaned];
-1 "Fitted NG monthly seasonal factors (Jan..Dec), as log-forward premium:";
show ([] month:`Jan`Feb`Mar`Apr`May`Jun`Jul`Aug`Sep`Oct`Nov`Dec; factor:fitted);
-1 "winter (Dec/Jan) avg=",(string avg fitted 11 0)," vs summer (Jun/Jul) avg=",string avg fitted 5 6;

/ Calibrate schwartz2 to an NG snapshot (mid-history; positive tenors only).
asofSnap:allDates (count allDates) div 3;
snap:.parser.futures.curveAt[longTable;asofSnap];
calCurve:select tenor,price from snap where tenor>0f, price>0f;
calRes:@[.commodity.calibrateCurve[;`schwartz2;()!()];calCurve;{[e] `status`errorMessage!(`ERROR;e)}];
-1 "";
$[`ERROR~calRes`status;
    -1 "NG schwartz2 snapshot calibration skipped: ",calRes`errorMessage;
    [-1 "NG schwartz2 fit as of ",(string asofSnap),": rmse=",string calRes`fitRmse;
     show calRes`calibratedParams]];
exit 0;
